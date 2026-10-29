"""Kalkulator laju transfer elektron non-adiabatik Marcus.

Modul ini mengubah parameter hasil DFT/MD menjadi estimasi laju transfer
elektron satu-elektron:

    k_ET = (2π/ħ) |H_AB|² / sqrt(4π λ k_B T)
           × exp[-(ΔG° + λ)² / (4 λ k_B T)]

Input dapat diberikan langsung (ΔG°, λ, H_AB) atau diturunkan dari:

- potensial reduksi donor dan akseptor (ΔG°),
- energi empat-titik DFT untuk reorganisasi inner-sphere (λ_i),
- model kontinum Marcus untuk reorganisasi outer-sphere (λ_o).

Catatan ilmiah penting: NAD+/NADP+ adalah akseptor **hidrida** dua-elektron
dalam metabolisme normal, bukan akseptor elektron luar-sphere satu-elektron.
Model Marcus di sini hanya sah untuk langkah transfer elektron satu-elektron
yang didefinisikan eksplisit. Untuk hidrida/PCET gunakan model yang sesuai.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JSON = PROJECT_ROOT / "results" / "marcus_electron_transfer.json"
DEFAULT_MD = PROJECT_ROOT / "results" / "marcus_electron_transfer.md"

KB_EV_PER_K = 8.617333262e-5  # eV/K
HBAR_EV_S = 6.582119569e-16  # eV·s
COULOMB_EV_ANGSTROM = 14.3996454784255  # e²/(4π ε₀) dalam eV·Å
CM1_PER_EV = 8065.544  # cm⁻¹ per eV
DEFAULT_TEMPERATURE_K = 298.15

# λ_outer = (e²/8π ε₀)(1/ε_op - 1/ε_s)(1/r_D + 1/r_A - 1/R)
OUTER_SPHERE_FACTOR = COULOMB_EV_ANGSTROM / 2.0


class MarcusError(RuntimeError):
    """Kesalahan input atau perhitungan Marcus."""


@dataclass(frozen=True)
class MarcusResult:
    delta_g_ev: float
    lambda_inner_ev: float
    lambda_outer_ev: float
    lambda_total_ev: float
    activation_energy_ev: float
    coupling_ev: float
    temperature_k: float
    prefactor_s_inv: float
    exponential_factor: float
    rate_s_inv: float

    def to_dict(self) -> dict[str, Any]:
        log10_rate = (
            math.log10(self.rate_s_inv) if self.rate_s_inv > 0 else None
        )
        return {
            "delta_g_ev": _round(self.delta_g_ev),
            "lambda_inner_ev": _round(self.lambda_inner_ev),
            "lambda_outer_ev": _round(self.lambda_outer_ev),
            "lambda_total_ev": _round(self.lambda_total_ev),
            "activation_energy_ev": _round(self.activation_energy_ev),
            "coupling_ev": _round(self.coupling_ev),
            "coupling_cm1": _round(self.coupling_ev * CM1_PER_EV),
            "temperature_k": _round(self.temperature_k),
            "thermal_energy_ev": _round(KB_EV_PER_K * self.temperature_k),
            "prefactor_s_inv": _scientific(self.prefactor_s_inv),
            "exponential_factor": _scientific(self.exponential_factor),
            "rate_s_inv": _scientific(self.rate_s_inv),
            "log10_rate_s_inv": _round(log10_rate)
            if log10_rate is not None
            else None,
        }


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _round(value: float, digits: int = 6) -> float:
    return round(float(value), digits)


def _scientific(value: float, digits: int = 6) -> float:
    """Bulatkan ke digit signifikan agar laju sangat kecil tidak jadi 0."""
    return float(f"{float(value):.{digits}e}")


def _read_json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise MarcusError(f"File tidak ditemukan: {path}") from exc
    except json.JSONDecodeError as exc:
        raise MarcusError(f"JSON tidak valid: {path} ({exc})") from exc
    if not isinstance(data, dict):
        raise MarcusError(f"JSON harus berupa objek: {path}")
    return data


def _positive(value: float, label: str) -> float:
    if not math.isfinite(value) or value <= 0:
        raise MarcusError(f"{label} harus > 0; diterima {value}.")
    return float(value)


def marcus_rate(
    delta_g_ev: float,
    lambda_total_ev: float,
    coupling_ev: float,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
) -> MarcusResult:
    """Hitung laju Marcus untuk satu donor-akseptor pada jarak tetap."""
    lambda_total_ev = _positive(lambda_total_ev, "lambda_total_ev")
    coupling_ev = float(coupling_ev)
    temperature_k = _positive(temperature_k, "temperature_k")
    if not math.isfinite(delta_g_ev):
        raise MarcusError("delta_g_ev harus berupa angka hingga.")

    thermal = KB_EV_PER_K * temperature_k
    activation = (delta_g_ev + lambda_total_ev) ** 2 / (
        4.0 * lambda_total_ev
    )
    prefactor = (
        (2.0 * math.pi / HBAR_EV_S)
        * coupling_ev**2
        / math.sqrt(4.0 * math.pi * lambda_total_ev * thermal)
    )
    exponential = math.exp(-activation / thermal)
    rate = prefactor * exponential
    if not math.isfinite(prefactor) or not math.isfinite(rate):
        raise MarcusError(
            "Hasil laju tidak hingga; periksa coupling dan lambda."
        )
    return MarcusResult(
        delta_g_ev=float(delta_g_ev),
        lambda_inner_ev=0.0,
        lambda_outer_ev=0.0,
        lambda_total_ev=lambda_total_ev,
        activation_energy_ev=activation,
        coupling_ev=coupling_ev,
        temperature_k=temperature_k,
        prefactor_s_inv=prefactor,
        exponential_factor=exponential,
        rate_s_inv=rate,
    )


def delta_g_from_reduction_potentials(
    donor_e0_v: float,
    acceptor_e0_v: float,
    n_electrons: int = 1,
) -> float:
    """ΔG° (eV) untuk D + A -> D+ + A− dari potensial reduksi.

    E°(D+/D) dan E°(A/A−) harus memakai referensi elektroda yang sama.
    Jika yang tersedia adalah potensial oksidasi E_ox(D), gunakan
    donor_e0_v = -E_ox(D).
    """
    if n_electrons < 1:
        raise MarcusError("n_electrons harus >= 1.")
    return float(n_electrons) * (float(donor_e0_v) - float(acceptor_e0_v))


def outer_sphere_lambda(
    donor_radius_angstrom: float,
    acceptor_radius_angstrom: float,
    distance_angstrom: float,
    dielectric_static: float,
    dielectric_optical: float,
) -> float:
    """Reorganisasi outer-sphere Marcus untuk dua bola dalam kontinum."""
    r_donor = _positive(donor_radius_angstrom, "donor_radius_angstrom")
    r_acceptor = _positive(
        acceptor_radius_angstrom, "acceptor_radius_angstrom"
    )
    distance = _positive(distance_angstrom, "distance_angstrom")
    eps_s = _positive(dielectric_static, "dielectric_static")
    eps_op = _positive(dielectric_optical, "dielectric_optical")
    if eps_s <= 1.0:
        raise MarcusError("dielectric_static harus > 1.")
    if eps_op >= eps_s:
        raise MarcusError(
            "dielectric_optical harus lebih kecil dari dielectric_static."
        )
    if distance < r_donor + r_acceptor - 1e-9:
        raise MarcusError(
            "distance_angstrom harus >= jumlah radius donor dan akseptor."
        )
    geometric = 1.0 / r_donor + 1.0 / r_acceptor - 1.0 / distance
    if geometric <= 0:
        raise MarcusError(
            "Geometri radius/jarak menghasilkan faktor reorganisasi <= 0."
        )
    value = (
        OUTER_SPHERE_FACTOR
        * (1.0 / eps_op - 1.0 / eps_s)
        * geometric
    )
    if value <= 0 or not math.isfinite(value):
        raise MarcusError("lambda outer-sphere tidak valid.")
    return value


def four_point_reorganization(data: dict[str, Any]) -> dict[str, Any]:
    """Hitung λ inner-sphere dari energi empat-titik DFT.

    Skema JSON::

        {
          "donor": {
            "neutral_at_neutral_eV": 0.0,
            "cation_at_neutral_eV": 7.00,
            "cation_at_cation_eV": 6.80,
            "neutral_at_cation_eV": 0.10
          },
          "acceptor": {
            "neutral_at_neutral_eV": 0.0,
            "anion_at_neutral_eV": -1.00,
            "anion_at_anion_eV": -1.25,
            "neutral_at_anion_eV": 0.10
          }
        }
    """
    try:
        donor = data["donor"]
        acceptor = data["acceptor"]
    except (KeyError, TypeError) as exc:
        raise MarcusError(
            "four_point harus memiliki objek 'donor' dan 'acceptor'."
        ) from exc

    def energy(block: dict[str, Any], key: str, label: str) -> float:
        if key not in block:
            raise MarcusError(f"four_point.{label}.{key} tidak ditemukan.")
        try:
            value = float(block[key])
        except (TypeError, ValueError) as exc:
            raise MarcusError(
                f"four_point.{label}.{key} bukan angka."
            ) from exc
        if not math.isfinite(value):
            raise MarcusError(f"four_point.{label}.{key} tidak hingga.")
        return value

    d_nn = energy(donor, "neutral_at_neutral_eV", "donor")
    d_cn = energy(donor, "cation_at_neutral_eV", "donor")
    d_cc = energy(donor, "cation_at_cation_eV", "donor")
    d_nc = energy(donor, "neutral_at_cation_eV", "donor")
    a_nn = energy(acceptor, "neutral_at_neutral_eV", "acceptor")
    a_an = energy(acceptor, "anion_at_neutral_eV", "acceptor")
    a_aa = energy(acceptor, "anion_at_anion_eV", "acceptor")
    a_na = energy(acceptor, "neutral_at_anion_eV", "acceptor")

    lambda_donor = (d_cn - d_cc) + (d_nc - d_nn)
    lambda_acceptor = (a_an - a_aa) + (a_na - a_nn)
    if lambda_donor < -1e-6 or lambda_acceptor < -1e-6:
        raise MarcusError(
            "lambda four-point negatif; periksa konsistensi level teori dan "
            "koreksi ZPE/solvasi."
        )
    lambda_donor = max(0.0, lambda_donor)
    lambda_acceptor = max(0.0, lambda_acceptor)
    return {
        "lambda_donor_ev": lambda_donor,
        "lambda_acceptor_ev": lambda_acceptor,
        "lambda_inner_ev": lambda_donor + lambda_acceptor,
        "components": {
            "donor_vertical_relaxation_ev": d_cn - d_cc,
            "donor_geometry_relaxation_ev": d_nc - d_nn,
            "acceptor_vertical_relaxation_ev": a_an - a_aa,
            "acceptor_geometry_relaxation_ev": a_na - a_nn,
        },
    }


def _coupling_to_ev(value: float, unit: str) -> float:
    normalized = unit.strip().lower()
    if normalized in {"ev", "electronvolt", "electron-volt"}:
        return float(value)
    if normalized in {"mev", "millielectronvolt"}:
        return float(value) / 1000.0
    if normalized in {"cm-1", "cm^-1", "wavenumber"}:
        return float(value) / CM1_PER_EV
    raise MarcusError(f"Satuan coupling tidak dikenal: {unit}")


def _get(
    args: argparse.Namespace,
    config: dict[str, Any],
    attribute: str,
    key: str,
    default: Any = None,
) -> Any:
    value = getattr(args, attribute, None)
    if value is not None:
        return value
    return config.get(key, default)


def _resolve_four_point(
    args: argparse.Namespace, config: dict[str, Any]
) -> tuple[dict[str, Any] | None, str | None]:
    if args.four_point_json:
        path = Path(args.four_point_json).expanduser().resolve()
        return _read_json(path), str(path)
    data = config.get("four_point")
    if isinstance(data, dict):
        return data, "config.four_point"
    return None, None


def _build_result(
    args: argparse.Namespace,
) -> tuple[MarcusResult, dict[str, Any], list[str]]:
    config = (
        _read_json(Path(args.config).expanduser().resolve())
        if args.config
        else {}
    )
    warnings: list[str] = []

    delta_g = _get(args, config, "delta_g_ev", "delta_g_ev")
    donor_e0 = _get(args, config, "donor_e0_v", "donor_e0_v")
    acceptor_e0 = _get(args, config, "acceptor_e0_v", "acceptor_e0_v")
    n_electrons_value = _get(
        args, config, "n_electrons", "n_electrons", 1
    )
    n_electrons = int(
        1 if n_electrons_value is None else n_electrons_value
    )
    if n_electrons < 1:
        raise MarcusError("n_electrons harus >= 1.")
    if delta_g is None:
        if donor_e0 is None or acceptor_e0 is None:
            raise MarcusError(
                "Berikan --delta-g-ev, atau --donor-e0-v dan --acceptor-e0-v."
            )
        delta_g = delta_g_from_reduction_potentials(
            float(donor_e0), float(acceptor_e0), n_electrons
        )
    elif donor_e0 is not None and acceptor_e0 is not None:
        check = delta_g_from_reduction_potentials(
            float(donor_e0), float(acceptor_e0), n_electrons
        )
        if abs(check - float(delta_g)) > 1e-6:
            warnings.append(
                "delta_g_ev berbeda dari nilai yang dihitung dari potensial; "
                "nilai delta_g_ev dipakai."
            )
    delta_g = float(delta_g)
    if n_electrons != 1:
        warnings.append(
            f"n_electrons={n_electrons} hanya dipakai untuk delta-G; persamaan "
            "Marcus yang dihitung tetap satu-elektron."
        )

    four_point, four_point_source = _resolve_four_point(args, config)
    lambda_inner_cli = _get(
        args, config, "lambda_inner_ev", "lambda_inner_ev"
    )
    lambda_inner: float
    four_point_payload: dict[str, Any] | None = None
    if four_point is not None:
        four_point_payload = four_point_reorganization(four_point)
        lambda_inner = float(four_point_payload["lambda_inner_ev"])
        if lambda_inner_cli is not None:
            warnings.append(
                "lambda_inner_ev dari CLI diabaikan karena four_point dipakai."
            )
    else:
        if lambda_inner_cli is None:
            raise MarcusError(
                "Berikan --lambda-inner-ev atau --four-point-json."
            )
        lambda_inner = float(lambda_inner_cli)
    if lambda_inner < 0:
        raise MarcusError("lambda_inner_ev tidak boleh negatif.")

    lambda_outer_cli = _get(
        args, config, "lambda_outer_ev", "lambda_outer_ev"
    )
    outer_params = {
        "donor_radius_angstrom": _get(
            args, config, "donor_radius_a", "donor_radius_angstrom"
        ),
        "acceptor_radius_angstrom": _get(
            args, config, "acceptor_radius_a", "acceptor_radius_angstrom"
        ),
        "distance_angstrom": _get(
            args, config, "distance_a", "distance_angstrom"
        ),
        "dielectric_static": _get(
            args, config, "eps_static", "dielectric_static"
        ),
        "dielectric_optical": _get(
            args, config, "eps_optical", "dielectric_optical"
        ),
    }
    if lambda_outer_cli is not None:
        lambda_outer = float(lambda_outer_cli)
        if all(value is not None for value in outer_params.values()):
            warnings.append(
                "lambda_outer_ev dari CLI dipakai; parameter bola/kontinum "
                "diabaikan."
            )
    else:
        missing = [
            label
            for label, value in outer_params.items()
            if value is None
        ]
        if missing:
            raise MarcusError(
                "Berikan --lambda-outer-ev, atau lengkapi parameter "
                "outer-sphere: " + ", ".join(missing)
            )
        lambda_outer = outer_sphere_lambda(
            float(outer_params["donor_radius_angstrom"]),
            float(outer_params["acceptor_radius_angstrom"]),
            float(outer_params["distance_angstrom"]),
            float(outer_params["dielectric_static"]),
            float(outer_params["dielectric_optical"]),
        )
    if lambda_outer < 0:
        raise MarcusError("lambda_outer_ev tidak boleh negatif.")

    lambda_total = lambda_inner + lambda_outer
    if lambda_total <= 0:
        raise MarcusError("lambda_total harus > 0.")

    coupling_raw = _get(args, config, "coupling", "coupling_ev")
    if coupling_raw is None:
        coupling_raw = config.get("coupling")
    if coupling_raw is None:
        raise MarcusError("Berikan --coupling (eV/meV/cm-1).")
    coupling_unit = _get(
        args, config, "coupling_unit", "coupling_unit", "eV"
    )
    coupling = _coupling_to_ev(float(coupling_raw), str(coupling_unit))
    if not math.isfinite(coupling):
        raise MarcusError("Coupling harus berupa angka hingga.")
    if coupling < 0:
        warnings.append("Coupling negatif; yang dipakai adalah |H_AB|^2.")
        coupling = abs(coupling)
    if coupling == 0:
        warnings.append("H_AB = 0; laju Marcus menjadi nol.")

    temperature_value = _get(
        args,
        config,
        "temperature_k",
        "temperature_k",
        DEFAULT_TEMPERATURE_K,
    )
    temperature = float(
        DEFAULT_TEMPERATURE_K
        if temperature_value is None
        else temperature_value
    )
    result = marcus_rate(delta_g, lambda_total, coupling, temperature)
    result = MarcusResult(
        delta_g_ev=result.delta_g_ev,
        lambda_inner_ev=lambda_inner,
        lambda_outer_ev=lambda_outer,
        lambda_total_ev=lambda_total,
        activation_energy_ev=result.activation_energy_ev,
        coupling_ev=result.coupling_ev,
        temperature_k=result.temperature_k,
        prefactor_s_inv=result.prefactor_s_inv,
        exponential_factor=result.exponential_factor,
        rate_s_inv=result.rate_s_inv,
    )

    inputs = {
        "delta_g_ev": _round(delta_g),
        "donor_e0_v": _round(float(donor_e0)) if donor_e0 is not None else None,
        "acceptor_e0_v": _round(float(acceptor_e0))
        if acceptor_e0 is not None
        else None,
        "n_electrons": n_electrons,
        "lambda_inner_ev": _round(lambda_inner),
        "lambda_outer_ev": _round(lambda_outer),
        "lambda_source": {
            "inner": four_point_source or "direct",
            "outer": "direct" if lambda_outer_cli is not None else "continuum",
        },
        "four_point": four_point_payload,
        "outer_sphere": {
            key: _round(float(value))
            for key, value in outer_params.items()
            if value is not None
        }
        or None,
        "coupling_ev": _round(coupling),
        "coupling_cm1": _round(coupling * CM1_PER_EV),
        "coupling_unit": coupling_unit,
        "temperature_k": _round(temperature),
        "thermal_energy_ev": _round(KB_EV_PER_K * temperature),
    }
    return result, inputs, warnings


def _write_markdown(
    path: Path,
    *,
    result: MarcusResult,
    inputs: dict[str, Any],
    warnings: list[str],
) -> None:
    rate = result.to_dict()
    log_rate = (
        f"{rate['log10_rate_s_inv']:.2f}"
        if rate["log10_rate_s_inv"] is not None
        else "tidak terdefinisi (k = 0)"
    )
    lines = [
        "# Estimasi Laju Transfer Elektron Marcus",
        "",
        f"> Dibuat: {_now()}",
        "> Model non-adiabatik Marcus untuk transfer elektron satu-elektron.",
        "",
        "## Parameter",
        "",
        "| Parameter | Nilai |",
        "|---|---|",
        f"| ΔG° | {inputs['delta_g_ev']:.6f} eV |",
        f"| λ inner-sphere | {inputs['lambda_inner_ev']:.6f} eV |",
        f"| λ outer-sphere | {inputs['lambda_outer_ev']:.6f} eV |",
        f"| λ total | {result.lambda_total_ev:.6f} eV |",
        f"| H_AB | {inputs['coupling_ev']:.6f} eV "
        f"({inputs['coupling_cm1']:.3f} cm⁻¹) |",
        f"| Temperatur | {inputs['temperature_k']:.2f} K |",
        f"| k_B T | {inputs['thermal_energy_ev']:.6f} eV |",
        "",
        "## Hasil",
        "",
        "| Besaran | Nilai |",
        "|---|---|",
        f"| Energi aktivasi ΔG‡ | {result.activation_energy_ev:.6f} eV |",
        f"| Prefaktor | {result.prefactor_s_inv:.6e} s⁻¹ |",
        f"| Faktor eksponensial | {result.exponential_factor:.6e} |",
        f"| Laju k_ET | {result.rate_s_inv:.6e} s⁻¹ |",
        f"| log10(k_ET / s⁻¹) | {log_rate} |",
        "",
        "Persamaan yang dipakai:",
        "",
        "```text",
        "k_ET = (2π/ħ) |H_AB|² / sqrt(4π λ k_B T)",
        "       × exp[-(ΔG° + λ)² / (4 λ k_B T)]",
        "```",
        "",
        "## Asumsi dan Keterbatasan",
        "",
        "- Model ini adalah Marcus non-adiabatik satu-elektron untuk pasangan "
        "donor-akseptor pada jarak dan orientasi tetap.",
        "- NAD+/NADP+ secara fisiologis menerima **hidrida** (2e⁻/H⁺) pada "
        "cincin nikotinamida; transfer elektron satu-elektron ke NAD+ "
        "menghasilkan radikal NAD• yang reaktif. Model ini tidak otomatis "
        "mewakili jalur hidrida.",
        "- Oksidasi DHN bersifat pH-dependent dan dapat berlangsung melalui "
        "PCET (deprotonasi + oksidasi), bukan ET luar-sphere sederhana.",
        "- λ inner-sphere dari four-point DFT bergantung pada fungsi "
        "fungsional, basis set, koreksi ZPE, dan model solvasi.",
        "- H_AB dari DFT sangat sensitif terhadap jarak, orientasi, dan "
        "tingkat teori; nilai kecil (≤10 meV) berada di batas akurasi.",
        "- Persamaan klasik mengabaikan nuclear tunneling, efek kuantum "
        "ikatan hidrogen, dan dinamika pelarut non-kontinum.",
        "- k_ET adalah laju unimolekuler untuk pasangan pada jarak tertentu, "
        "bukan konstanta bimolekuler; untuk larutan encer perlu koreksi "
        "difusi/encounter.",
        "- Angka ini adalah estimasi komputasi, bukan bukti eksperimental.",
        "",
    ]
    if warnings:
        lines.extend(["## Peringatan", ""])
        lines.extend(f"- {warning}" for warning in warnings)
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def _self_test() -> None:
    result = marcus_rate(-0.5, 1.0, 0.001, 298.15)
    expected_activation = 0.5**2 / (4.0 * 1.0)
    if not math.isclose(
        result.activation_energy_ev, expected_activation, rel_tol=1e-12
    ):
        raise MarcusError("Self-test dG# gagal.")
    expected_rate = (
        (2.0 * math.pi / HBAR_EV_S)
        * 0.001**2
        / math.sqrt(4.0 * math.pi * 1.0 * KB_EV_PER_K * 298.15)
        * math.exp(-expected_activation / (KB_EV_PER_K * 298.15))
    )
    if not math.isclose(result.rate_s_inv, expected_rate, rel_tol=1e-12):
        raise MarcusError("Self-test laju gagal.")
    outer = outer_sphere_lambda(3.0, 3.0, 6.0, 78.355, 1.776)
    expected_outer = (
        OUTER_SPHERE_FACTOR
        * (1.0 / 1.776 - 1.0 / 78.355)
        * (1.0 / 3.0 + 1.0 / 3.0 - 1.0 / 6.0)
    )
    if not math.isclose(outer, expected_outer, rel_tol=1e-12):
        raise MarcusError("Self-test outer-sphere gagal.")
    four_point = four_point_reorganization(
        {
            "donor": {
                "neutral_at_neutral_eV": 0.0,
                "cation_at_neutral_eV": 7.0,
                "cation_at_cation_eV": 6.8,
                "neutral_at_cation_eV": 0.1,
            },
            "acceptor": {
                "neutral_at_neutral_eV": 0.0,
                "anion_at_neutral_eV": -1.0,
                "anion_at_anion_eV": -1.25,
                "neutral_at_anion_eV": 0.1,
            },
        }
    )
    if not math.isclose(
        four_point["lambda_inner_ev"], 0.65, rel_tol=1e-12
    ):
        raise MarcusError("Self-test four-point gagal.")
    if not math.isclose(
        delta_g_from_reduction_potentials(-0.5, -0.32), -0.18, rel_tol=1e-12
    ):
        raise MarcusError("Self-test potensial gagal.")
    print("[SELF-TEST] Semua pemeriksaan numerik lulus.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Estimasi laju transfer elektron Marcus dari parameter DFT/MD."
        )
    )
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--delta-g-ev", type=float, default=None)
    parser.add_argument("--donor-e0-v", type=float, default=None)
    parser.add_argument("--acceptor-e0-v", type=float, default=None)
    parser.add_argument("--n-electrons", type=int, default=None)
    parser.add_argument("--lambda-inner-ev", type=float, default=None)
    parser.add_argument("--lambda-outer-ev", type=float, default=None)
    parser.add_argument("--four-point-json", type=Path, default=None)
    parser.add_argument("--donor-radius-a", type=float, default=None)
    parser.add_argument("--acceptor-radius-a", type=float, default=None)
    parser.add_argument("--distance-a", type=float, default=None)
    parser.add_argument("--eps-static", type=float, default=None)
    parser.add_argument("--eps-optical", type=float, default=None)
    parser.add_argument("--coupling", type=float, default=None)
    parser.add_argument(
        "--coupling-unit",
        choices=["eV", "meV", "cm-1"],
        default=None,
        help="Satuan --coupling (default: eV; config: coupling_unit).",
    )
    parser.add_argument("--temperature-k", type=float, default=None)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--md", type=Path, default=DEFAULT_MD)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--example", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    # Hindari crash UnicodeEncodeError saat console Windows memakai cp1252.
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(errors="replace")
        except Exception:
            pass

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.self_test:
        try:
            _self_test()
            return 0
        except MarcusError as exc:
            print(f"[ERROR] {exc}")
            return 1

    if args.example:
        print(
            "Contoh:\n"
            "python src/marcus_electron_transfer.py \\\n"
            "  --delta-g-ev -0.18 \\\n"
            "  --lambda-inner-ev 0.45 --lambda-outer-ev 0.85 \\\n"
            "  --coupling 1.0 --coupling-unit meV"
        )
        return 0

    try:
        result, inputs, warnings = _build_result(args)
        payload = {
            "tool": "melanin_biotech_marcus_electron_transfer",
            "generated_at": _now(),
            "model": "non-adiabatic Marcus (one-electron, fixed distance)",
            "inputs": inputs,
            "results": result.to_dict(),
            "warnings": warnings,
            "caveats": [
                "NAD+/NADP+ is a hydride acceptor in normal metabolism; "
                "one-electron Marcus ET is a model, not the physiological "
                "hydride pathway.",
                "DHN oxidation is pH-dependent and may proceed by PCET.",
                "H_AB and lambda are highly sensitive to DFT functional, "
                "basis set, solvation, and geometry.",
                "k_ET is a unimolecular rate at a fixed donor-acceptor "
                "distance, not a bimolecular diffusion-corrected rate.",
            ],
        }
        json_path = Path(args.json).expanduser().resolve()
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        md_path = Path(args.md).expanduser().resolve()
        _write_markdown(md_path, result=result, inputs=inputs, warnings=warnings)

        print("[SUCCESS] Marcus electron transfer estimate")
        print(
            f"    DeltaG = {result.delta_g_ev:.4f} eV; "
            f"lambda = {result.lambda_total_ev:.4f} eV "
            f"(inner {result.lambda_inner_ev:.4f}, "
            f"outer {result.lambda_outer_ev:.4f})"
        )
        print(
            f"    DeltaG# = {result.activation_energy_ev:.4f} eV; "
            f"H_AB = {result.coupling_ev:.6f} eV"
        )
        print(f"    k_ET = {result.rate_s_inv:.4e} s^-1")
        if result.rate_s_inv > 0:
            print(f"    log10(k_ET) = {math.log10(result.rate_s_inv):.3f}")
        print(f"    JSON : {json_path}")
        print(f"    MD   : {md_path}")
        for warning in warnings:
            print(f"    [WARN] {warning}")
        return 0
    except MarcusError as exc:
        print(f"[ERROR] {exc}")
        return 2
    except Exception as exc:
        print(f"[ERROR] {type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
