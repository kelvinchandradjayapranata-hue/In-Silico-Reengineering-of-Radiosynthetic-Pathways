"""Model FBA radiosintesis melanin (versi 2).

Model ini adalah model stoikiometri tereduksi (bukan genome-scale) yang
menghubungkan input radiasi, metabolisme karbon pusat, dan dua jalur melanin:

1. Jalur DHN-melanin (NADPH/ATP-dependent):
   acetyl-CoA + 4 malonyl-CoA -> THN
   THN -> scytalone -> 1,3,8-THN -> vermelone -> DHN -> DHN-melanin

2. Jalur DOPA-melanin (tyrosinase, O2-dependent):
   L-tyrosine -> L-DOPA -> dopaquinone -> DOPA-melanin

Input radiasi dimodelkan sebagai:
- elektron radiasi -> NADH (RAD_NADH_SYN)
- energi radiasi -> ATP (RAD_ATP_SYN)

Caveat: koenzim (CoA, H+) dan reaksi transport tidak dimodelkan secara rinci;
model ini adalah kerangka stoikiometri konseptual, bukan prediksi kuantitatif
yang tervalidasi eksperimental.
"""

from __future__ import annotations

import argparse
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from cobra import Metabolite, Model, Reaction
from cobra.io import write_sbml_model

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / "results"
MODELS_DIR = PROJECT_ROOT / "models"
FBA_RESULTS_PATH = RESULTS_DIR / "metabolic_fba.json"
SENSITIVITY_RESULTS_PATH = RESULTS_DIR / "fba_sensitivity.json"
SBML_PATH = MODELS_DIR / "melanin_radiosynthesis_fba_v2.xml"

MODEL_ID = "melanin_radiosynthesis_fba_v2"

# Fluks yang ditampilkan di JSON `fluxes` dan grafik. Semua reaksi tetap
# tersedia di `all_fluxes`.
KEY_FLUX_IDS: tuple[str, ...] = (
    "EX_GLC",
    "GLYCOLYSIS",
    "ACCARBOXYLASE",
    "PKS_THN",
    "MELANIN_DHN_EXPORT",
    "TYROSINASE_CRESOLASE",
    "MELANIN_DOPA_EXPORT",
    "RESPIRATION",
    "TRANS_HYDROGENASE",
    "RAD_ELECTRON_IN",
    "RAD_ENERGY_IN",
    "ATPM",
    "BIOMASS",
)

ASSUMPTIONS: tuple[str, ...] = (
    "Model adalah stoikiometri tereduksi, bukan genome-scale model.",
    "1 acetyl-CoA + 4 malonyl-CoA menghasilkan 1 THN (PKS non-reducing, lumped).",
    "Jalur DHN: 2 NADPH per DHN (THN reductase dan 1,3,8-THN reductase).",
    "1 ATP dan 1 CO2 dikonsumsi per malonyl-CoA (acetyl-CoA carboxylase).",
    "1 glucose -> 2 acetyl-CoA + 2 NADH + 2 ATP + 2 CO2 (glycolysis lumped).",
    "1 glucose -> 1 pentose phosphate + 2 NADPH + 1 CO2 (PPP oksidatif lumped).",
    "Respirasi lumped: 2 NADH + 0.5 O2 -> 2 NAD+ + 2 ATP (P/O = 1).",
    "Radiasi dimodelkan sebagai input elektron (-> NADH) dan energi (-> ATP), "
    "bukan model fotokimia/kuantum.",
    "Koenzim seperti CoA dan H+ dihilangkan; neraca massa/charge tidak lengkap.",
    "Reaksi biomassa adalah placeholder, bukan komposisi biomassa terukur.",
    "Jalur DOPA adalah reaksi tyrosinase + polimerisasi yang dilumpukan.",
    "Semua hasil adalah estimasi in silico, bukan validasi eksperimental.",
)


def _metabolite(metabolite_id: str, name: str) -> Metabolite:
    return Metabolite(metabolite_id, name=name, compartment="c")


def _add_reaction(
    model: Model,
    reaction_id: str,
    name: str,
    stoichiometry: dict[Metabolite, float],
    lower_bound: float = 0.0,
    upper_bound: float = 1000.0,
) -> Reaction:
    reaction = Reaction(reaction_id)
    reaction.name = name
    reaction.lower_bound = lower_bound
    reaction.upper_bound = upper_bound
    reaction.add_metabolites(stoichiometry)
    model.add_reactions([reaction])
    return reaction


def build_radiosynthesis_model(
    electron_input_limit: float = 10.0,
    energy_input_limit: float = 10.0,
    glucose_limit: float = 5.0,
    oxygen_limit: float = 20.0,
    tyrosine_limit: float = 1.0,
    atp_maintenance: float = 10.0,
    biomass_requirement: float = 0.1,
) -> Model:
    """Bangun model FBA radiosintesis melanin versi 2.

    Parameter input adalah batas uptake (mmol gDW^-1 h^-1) dan kebutuhan
    maintenance/biomassa. Nilai default dipilih agar energi radiasi menjadi
    bagian dari neraca ATP, tetapi tidak memaksa semua input radiasi terpakai.
    """
    for name, value in {
        "electron_input_limit": electron_input_limit,
        "energy_input_limit": energy_input_limit,
        "glucose_limit": glucose_limit,
        "oxygen_limit": oxygen_limit,
        "tyrosine_limit": tyrosine_limit,
        "atp_maintenance": atp_maintenance,
        "biomass_requirement": biomass_requirement,
    }.items():
        if value < 0:
            raise ValueError(f"{name} tidak boleh negatif: {value}")

    model = Model(MODEL_ID)
    model.name = "Melanin radiosynthesis reduced stoichiometric model v2"

    # --- Metabolisme karbon pusat ---
    glucose = _metabolite("glc_c", "D-glucose")
    oxygen = _metabolite("o2_c", "Oxygen")
    tyrosine = _metabolite("tyr_c", "L-tyrosine")
    co2 = _metabolite("co2_c", "CO2")
    water = _metabolite("h2o_c", "H2O")
    acetyl_coa = _metabolite("accoa_c", "Acetyl-CoA")
    malonyl_coa = _metabolite("malcoa_c", "Malonyl-CoA")
    pentose_phosphate = _metabolite("ppp_c", "Pentose phosphate (lumped)")

    # --- Redoks dan energi ---
    nad = _metabolite("nad_c", "NAD+")
    nadh = _metabolite("nadh_c", "NADH")
    nadp = _metabolite("nadp_c", "NADP+")
    nadph = _metabolite("nadph_c", "NADPH")
    atp = _metabolite("atp_c", "ATP")
    adp = _metabolite("adp_c", "ADP")
    phosphate = _metabolite("pi_c", "Phosphate")
    radiation_electron = _metabolite("e_red_c", "Radiation-derived electron")
    radiation_energy = _metabolite("e_energy_c", "Radiation-derived energy")

    # --- Jalur DHN-melanin ---
    thn = _metabolite("thn_c", "1,3,6,8-tetrahydroxynaphthalene (THN)")
    scytalone = _metabolite("scytalone_c", "Scytalone")
    trihydroxynaphthalene = _metabolite(
        "t3hn_c", "1,3,8-trihydroxynaphthalene (1,3,8-THN)"
    )
    vermelone = _metabolite("vermelone_c", "Vermelone")
    dhn = _metabolite("dhn_c", "1,8-dihydroxynaphthalene (DHN)")
    dhn_melanin = _metabolite("melanin_dhn_c", "DHN-melanin")

    # --- Jalur DOPA-melanin ---
    dopa = _metabolite("dopa_c", "L-DOPA")
    dopaquinone = _metabolite("dopaquinone_c", "Dopaquinone")
    dopa_melanin = _metabolite("melanin_dopa_c", "DOPA-melanin")

    # --- Biomassa placeholder ---
    biomass = _metabolite("biomass_c", "Biomass precursor pool (placeholder)")

    # --- Uptake dari medium ---
    _add_reaction(
        model,
        "EX_GLC",
        "Glucose uptake",
        {glucose: 1.0},
        upper_bound=glucose_limit,
    )
    _add_reaction(
        model,
        "EX_O2",
        "Oxygen uptake",
        {oxygen: 1.0},
        upper_bound=oxygen_limit,
    )
    _add_reaction(
        model,
        "EX_TYR",
        "L-tyrosine uptake",
        {tyrosine: 1.0},
        upper_bound=tyrosine_limit,
    )

    # --- Metabolisme karbon pusat ---
    _add_reaction(
        model,
        "GLYCOLYSIS",
        "Glycolysis (lumped): glucose -> 2 acetyl-CoA + 2 NADH + 2 ATP",
        {
            glucose: -1.0,
            nad: -2.0,
            adp: -2.0,
            phosphate: -2.0,
            acetyl_coa: 2.0,
            nadh: 2.0,
            atp: 2.0,
            co2: 2.0,
            water: 2.0,
        },
    )
    _add_reaction(
        model,
        "PPP",
        "Pentose phosphate pathway (lumped NADPH source)",
        {
            glucose: -1.0,
            nadp: -2.0,
            water: -1.0,
            pentose_phosphate: 1.0,
            nadph: 2.0,
            co2: 1.0,
        },
    )
    _add_reaction(
        model,
        "PPP_EXPORT",
        "Pentose phosphate overflow",
        {pentose_phosphate: -1.0},
    )
    _add_reaction(
        model,
        "ACCARBOXYLASE",
        "Acetyl-CoA carboxylase: acetyl-CoA + ATP + CO2 -> malonyl-CoA",
        {
            acetyl_coa: -1.0,
            atp: -1.0,
            co2: -1.0,
            malonyl_coa: 1.0,
            adp: 1.0,
            phosphate: 1.0,
        },
    )

    # --- Jalur DHN-melanin ---
    _add_reaction(
        model,
        "PKS_THN",
        "Non-reducing PKS: acetyl-CoA + 4 malonyl-CoA -> THN",
        {
            acetyl_coa: -1.0,
            malonyl_coa: -4.0,
            thn: 1.0,
            co2: 4.0,
        },
    )
    _add_reaction(
        model,
        "THN_REDUCTASE",
        "THN reductase: THN + NADPH -> scytalone",
        {thn: -1.0, nadph: -1.0, scytalone: 1.0, nadp: 1.0},
    )
    _add_reaction(
        model,
        "SCYTALONE_DEHYDRATASE",
        "Scytalone dehydratase: scytalone -> 1,3,8-THN",
        {scytalone: -1.0, trihydroxynaphthalene: 1.0, water: 1.0},
    )
    _add_reaction(
        model,
        "T3HN_REDUCTASE",
        "1,3,8-THN reductase: 1,3,8-THN + NADPH -> vermelone",
        {
            trihydroxynaphthalene: -1.0,
            nadph: -1.0,
            vermelone: 1.0,
            nadp: 1.0,
        },
    )
    _add_reaction(
        model,
        "VERMELONE_DEHYDRATASE",
        "Vermelone dehydratase: vermelone -> DHN",
        {vermelone: -1.0, dhn: 1.0, water: 1.0},
    )
    _add_reaction(
        model,
        "LACCASE_DHN",
        "Laccase oxidation/polymerization: DHN -> DHN-melanin",
        {dhn: -1.0, oxygen: -0.5, dhn_melanin: 1.0, water: 1.0},
    )
    dhn_export = _add_reaction(
        model,
        "MELANIN_DHN_EXPORT",
        "DHN-melanin export",
        {dhn_melanin: -1.0},
    )

    # --- Jalur DOPA-melanin ---
    _add_reaction(
        model,
        "TYROSINASE_CRESOLASE",
        "Tyrosinase cresolase: L-tyrosine + O2 -> L-DOPA",
        {tyrosine: -1.0, oxygen: -1.0, dopa: 1.0, water: 1.0},
    )
    _add_reaction(
        model,
        "TYROSINASE_CATECHOLASE",
        "Tyrosinase catecholase: L-DOPA + 0.5 O2 -> dopaquinone",
        {dopa: -1.0, oxygen: -0.5, dopaquinone: 1.0, water: 1.0},
    )
    _add_reaction(
        model,
        "DOPA_MELANIN_FORMATION",
        "Dopaquinone polymerization -> DOPA-melanin (lumped)",
        {dopaquinone: -1.0, dopa_melanin: 1.0},
    )
    dopa_export = _add_reaction(
        model,
        "MELANIN_DOPA_EXPORT",
        "DOPA-melanin export",
        {dopa_melanin: -1.0},
    )

    # --- Respirasi, transhidrogenase, maintenance ---
    _add_reaction(
        model,
        "RESPIRATION",
        "Oxidative phosphorylation (lumped, P/O = 1): 2 NADH + 0.5 O2 -> 2 ATP",
        {
            nadh: -2.0,
            oxygen: -0.5,
            adp: -2.0,
            phosphate: -2.0,
            nad: 2.0,
            atp: 2.0,
            water: 2.0,
        },
    )
    _add_reaction(
        model,
        "TRANS_HYDROGENASE",
        "NADH/NADPH transhydrogenase",
        {nadh: -1.0, nadp: -1.0, nad: 1.0, nadph: 1.0},
        lower_bound=-1000.0,
        upper_bound=1000.0,
    )
    _add_reaction(
        model,
        "ATPM",
        "ATP maintenance requirement",
        {atp: -1.0, water: -1.0, adp: 1.0, phosphate: 1.0},
        lower_bound=atp_maintenance,
        upper_bound=atp_maintenance,
    )
    _add_reaction(
        model,
        "BIOMASS",
        "Biomass precursor demand (placeholder)",
        {
            acetyl_coa: -2.0,
            pentose_phosphate: -1.0,
            atp: -5.0,
            nadph: -2.0,
            biomass: 1.0,
            adp: 5.0,
            phosphate: 5.0,
            nadp: 2.0,
        },
        lower_bound=biomass_requirement,
        upper_bound=biomass_requirement,
    )
    _add_reaction(
        model,
        "BIOMASS_EXPORT",
        "Biomass sink",
        {biomass: -1.0},
    )

    # --- Input radiasi ---
    _add_reaction(
        model,
        "RAD_ELECTRON_IN",
        "Radiation electron input",
        {radiation_electron: 1.0},
        upper_bound=electron_input_limit,
    )
    _add_reaction(
        model,
        "RAD_ENERGY_IN",
        "Radiation energy input",
        {radiation_energy: 1.0},
        upper_bound=energy_input_limit,
    )
    _add_reaction(
        model,
        "RAD_NADH_SYN",
        "Radiation-driven NADH synthesis: electron + NAD+ -> NADH",
        {radiation_electron: -1.0, nad: -1.0, nadh: 1.0},
    )
    _add_reaction(
        model,
        "RAD_ATP_SYN",
        "Radiation-driven ATP synthesis: energy + ADP + Pi -> ATP",
        {
            radiation_energy: -1.0,
            adp: -1.0,
            phosphate: -1.0,
            atp: 1.0,
            water: 1.0,
        },
    )

    # --- Pertukaran CO2 dan H2O (reversible, agar neraca tidak mengikat) ---
    _add_reaction(
        model,
        "CO2_EX",
        "CO2 exchange",
        {co2: -1.0},
        lower_bound=-1000.0,
        upper_bound=1000.0,
    )
    _add_reaction(
        model,
        "H2O_EX",
        "H2O exchange",
        {water: -1.0},
        lower_bound=-1000.0,
        upper_bound=1000.0,
    )

    # Objective: total ekspor melanin (DHN + DOPA), bobot sama.
    model.objective = {dhn_export: 1.0, dopa_export: 1.0}
    model.objective_direction = "max"
    return model


def _safe_float(value: Any, digits: int = 6) -> float:
    number = float(value)
    if abs(number) < 1e-9:
        return 0.0
    return round(number, digits)


def _binding_constraints(fluxes: dict[str, float], model: Model) -> list[str]:
    constraints: list[str] = []
    checks = {
        "EX_GLC": "Glucose uptake",
        "EX_O2": "Oxygen uptake",
        "EX_TYR": "Tyrosine uptake",
        "RAD_ELECTRON_IN": "Radiation electron input",
        "RAD_ENERGY_IN": "Radiation energy input",
    }
    for reaction_id, label in checks.items():
        reaction = model.reactions.get_by_id(reaction_id)
        flux = fluxes.get(reaction_id, 0.0)
        if reaction.upper_bound > 0 and abs(flux - reaction.upper_bound) < 1e-6:
            constraints.append(f"{label} mencapai batas atas ({reaction.upper_bound:g})")
        if reaction.lower_bound < 0 and abs(flux - reaction.lower_bound) < 1e-6:
            constraints.append(f"{label} mencapai batas bawah ({reaction.lower_bound:g})")
    return constraints


def _warnings(fluxes: dict[str, float], inputs: dict[str, float]) -> list[str]:
    warnings: list[str] = []
    if fluxes["MELANIN_DHN_EXPORT"] <= 1e-9:
        warnings.append("Ekspor DHN-melanin nol pada kondisi ini.")
    if fluxes["MELANIN_DOPA_EXPORT"] <= 1e-9:
        warnings.append("Ekspor DOPA-melanin nol pada kondisi ini.")
    if fluxes["RAD_ELECTRON_IN"] <= 1e-9:
        warnings.append(
            "Input elektron radiasi tidak terpakai; NADH dari glikolisis dan PPP "
            "cukup untuk kebutuhan redoks."
        )
    if fluxes["RAD_ENERGY_IN"] <= 1e-9:
        warnings.append(
            "Input energi radiasi tidak terpakai; ATP dari glikolisis dan "
            "respirasi memenuhi kebutuhan."
        )
    if fluxes["RAD_ENERGY_IN"] >= inputs["radiation_energy_limit"] - 1e-6:
        warnings.append("Input energi radiasi menjadi faktor pembatas (binding).")
    if fluxes["RAD_ELECTRON_IN"] >= inputs["radiation_electron_limit"] - 1e-6:
        warnings.append("Input elektron radiasi menjadi faktor pembatas (binding).")
    return warnings


def run_fba(
    output_path: Path | None = FBA_RESULTS_PATH,
    sbml_path: Path | None = SBML_PATH,
    electron_input_limit: float = 10.0,
    energy_input_limit: float = 10.0,
    glucose_limit: float = 5.0,
    oxygen_limit: float = 20.0,
    tyrosine_limit: float = 1.0,
    atp_maintenance: float = 10.0,
    biomass_requirement: float = 0.1,
) -> dict[str, Any]:
    """Jalankan FBA, tulis JSON, dan (opsional) SBML."""
    inputs = {
        "glucose_uptake_limit": glucose_limit,
        "oxygen_uptake_limit": oxygen_limit,
        "tyrosine_uptake_limit": tyrosine_limit,
        "radiation_electron_limit": electron_input_limit,
        "radiation_energy_limit": energy_input_limit,
        "atp_maintenance": atp_maintenance,
        "biomass_requirement": biomass_requirement,
    }

    model = build_radiosynthesis_model(
        electron_input_limit=electron_input_limit,
        energy_input_limit=energy_input_limit,
        glucose_limit=glucose_limit,
        oxygen_limit=oxygen_limit,
        tyrosine_limit=tyrosine_limit,
        atp_maintenance=atp_maintenance,
        biomass_requirement=biomass_requirement,
    )
    solution = model.optimize()

    if solution.status != "optimal":
        raise RuntimeError(
            f"FBA tidak optimal (status: {solution.status}). Periksa batas input "
            "dan kebutuhan ATP/biomassa."
        )

    all_fluxes = {
        reaction.id: _safe_float(solution.fluxes[reaction.id])
        for reaction in model.reactions
    }
    key_fluxes = {reaction_id: all_fluxes[reaction_id] for reaction_id in KEY_FLUX_IDS}

    dhn_export = all_fluxes["MELANIN_DHN_EXPORT"]
    dopa_export = all_fluxes["MELANIN_DOPA_EXPORT"]
    glucose_used = all_fluxes["EX_GLC"]
    tyrosine_used = all_fluxes["EX_TYR"]
    electron_used = all_fluxes["RAD_ELECTRON_IN"]
    energy_used = all_fluxes["RAD_ENERGY_IN"]
    radiation_atp = all_fluxes["RAD_ATP_SYN"]

    atp_demand = (
        all_fluxes["ATPM"]
        + all_fluxes["ACCARBOXYLASE"]
        + 5.0 * all_fluxes["BIOMASS"]
    )
    radiation_atp_share = radiation_atp / atp_demand if atp_demand > 0 else 0.0

    result: dict[str, Any] = {
        "model": model.id,
        "model_version": "2.0",
        "status": solution.status,
        "objective": {
            "name": "total_melanin_export",
            "weights": {"MELANIN_DHN_EXPORT": 1.0, "MELANIN_DOPA_EXPORT": 1.0},
            "value": _safe_float(solution.objective_value),
        },
        "objective_value": _safe_float(solution.objective_value),
        "inputs": inputs,
        "fluxes": key_fluxes,
        "all_fluxes": all_fluxes,
        "pathway_fluxes": {
            "dhn_melanin_export": dhn_export,
            "dopa_melanin_export": dopa_export,
            "total_melanin_export": _safe_float(dhn_export + dopa_export),
        },
        "yields": {
            "dhn_melanin_per_glucose": _safe_float(
                dhn_export / glucose_used if glucose_used > 0 else 0.0
            ),
            "dhn_melanin_per_radiation_electron": _safe_float(
                dhn_export / electron_used if electron_used > 0 else 0.0
            ),
            "dhn_melanin_per_radiation_energy": _safe_float(
                dhn_export / energy_used if energy_used > 0 else 0.0
            ),
            "dopa_melanin_per_tyrosine": _safe_float(
                dopa_export / tyrosine_used if tyrosine_used > 0 else 0.0
            ),
        },
        "balance": {
            "glucose_used": glucose_used,
            "oxygen_used": all_fluxes["EX_O2"],
            "tyrosine_used": tyrosine_used,
            "co2_exported": all_fluxes["CO2_EX"],
            "radiation_electron_used": electron_used,
            "radiation_energy_used": energy_used,
            "radiation_atp_produced": radiation_atp,
            "radiation_atp_share_of_demand": _safe_float(radiation_atp_share),
        },
        "binding_constraints": _binding_constraints(all_fluxes, model),
        "warnings": _warnings(all_fluxes, inputs),
        "assumptions": list(ASSUMPTIONS),
        "interpretation": (
            "Model stoikiometri tereduksi untuk eksplorasi kualitatif konversi "
            "energi radiasi menjadi ATP/NADPH dan biosintesis melanin."
        ),
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }

    if output_path is not None:
        destination = Path(output_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temp_path = destination.with_suffix(destination.suffix + ".tmp")
        temp_path.write_text(
            json.dumps(result, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        temp_path.replace(destination)

    if sbml_path is not None:
        destination = Path(sbml_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        write_sbml_model(model, str(destination))
        result["sbml_path"] = str(destination.relative_to(PROJECT_ROOT))

    return result


def run_sensitivity(
    output_path: Path | None = SENSITIVITY_RESULTS_PATH,
    electron_limits: Sequence[float] = (0.0, 2.5, 5.0, 7.5, 10.0),
    energy_limits: Sequence[float] = (0.0, 2.5, 5.0, 7.5, 10.0),
    **model_kwargs: Any,
) -> dict[str, Any]:
    """Sweep batas input radiasi untuk melihat kondisi binding."""
    rows: list[dict[str, Any]] = []
    for electron_limit in electron_limits:
        for energy_limit in energy_limits:
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "ignore", message="Solver status is 'infeasible'"
                    )
                    result = run_fba(
                        output_path=None,
                        sbml_path=None,
                        electron_input_limit=electron_limit,
                        energy_input_limit=energy_limit,
                        **model_kwargs,
                    )
                rows.append(
                    {
                        "electron_limit": electron_limit,
                        "energy_limit": energy_limit,
                        "status": result["status"],
                        "objective_value": result["objective_value"],
                        "dhn_melanin_export": result["pathway_fluxes"][
                            "dhn_melanin_export"
                        ],
                        "dopa_melanin_export": result["pathway_fluxes"][
                            "dopa_melanin_export"
                        ],
                        "radiation_electron_used": result["balance"][
                            "radiation_electron_used"
                        ],
                        "radiation_energy_used": result["balance"][
                            "radiation_energy_used"
                        ],
                    }
                )
            except Exception as exc:
                rows.append(
                    {
                        "electron_limit": electron_limit,
                        "energy_limit": energy_limit,
                        "status": "error",
                        "error": str(exc),
                    }
                )

    payload: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "electron_limits": list(electron_limits),
        "energy_limits": list(energy_limits),
        "model_parameters": {
            "glucose_limit": model_kwargs.get("glucose_limit", 5.0),
            "oxygen_limit": model_kwargs.get("oxygen_limit", 20.0),
            "tyrosine_limit": model_kwargs.get("tyrosine_limit", 1.0),
            "atp_maintenance": model_kwargs.get("atp_maintenance", 10.0),
            "biomass_requirement": model_kwargs.get("biomass_requirement", 0.1),
        },
        "rows": rows,
        "interpretation": (
            "Sweep ini menunjukkan kondisi ketika input elektron/energi radiasi "
            "menjadi faktor pembatas. Nilai adalah estimasi in silico."
        ),
    }

    if output_path is not None:
        destination = Path(output_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temp_path = destination.with_suffix(destination.suffix + ".tmp")
        temp_path.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        temp_path.replace(destination)

    return payload


def _print_summary(result: dict[str, Any]) -> None:
    print("=== Hasil Simulasi FBA v2 ===")
    print(f"Status                : {result['status']}")
    print(f"Objective (total)     : {result['objective_value']}")
    print(f"DHN-melanin export    : {result['pathway_fluxes']['dhn_melanin_export']}")
    print(f"DOPA-melanin export   : {result['pathway_fluxes']['dopa_melanin_export']}")
    print(f"Glucose used          : {result['balance']['glucose_used']}")
    print(f"Radiation e- used     : {result['balance']['radiation_electron_used']}")
    print(f"Radiation energy used : {result['balance']['radiation_energy_used']}")
    print(
        "Radiation ATP share   : "
        f"{result['balance']['radiation_atp_share_of_demand']:.3%}"
    )
    if result["binding_constraints"]:
        print("Binding constraints   :")
        for item in result["binding_constraints"]:
            print(f"  - {item}")
    if result["warnings"]:
        print("Warnings              :")
        for item in result["warnings"]:
            print(f"  - {item}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Jalankan FBA radiosintesis melanin versi 2."
    )
    parser.add_argument("--electron-limit", type=float, default=10.0)
    parser.add_argument("--energy-limit", type=float, default=10.0)
    parser.add_argument("--glucose-limit", type=float, default=5.0)
    parser.add_argument("--oxygen-limit", type=float, default=20.0)
    parser.add_argument("--tyrosine-limit", type=float, default=1.0)
    parser.add_argument("--atp-maintenance", type=float, default=10.0)
    parser.add_argument("--biomass", type=float, default=0.1)
    parser.add_argument("--output", type=Path, default=FBA_RESULTS_PATH)
    parser.add_argument("--sbml", type=Path, default=SBML_PATH)
    parser.add_argument("--no-sbml", action="store_true")
    parser.add_argument(
        "--sensitivity",
        action="store_true",
        help="Jalankan sweep batas elektron/energi radiasi.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = run_fba(
        output_path=args.output,
        sbml_path=None if args.no_sbml else args.sbml,
        electron_input_limit=args.electron_limit,
        energy_input_limit=args.energy_limit,
        glucose_limit=args.glucose_limit,
        oxygen_limit=args.oxygen_limit,
        tyrosine_limit=args.tyrosine_limit,
        atp_maintenance=args.atp_maintenance,
        biomass_requirement=args.biomass,
    )
    _print_summary(result)
    print(f"\n[SUCCESS] Hasil FBA disimpan ke: {args.output}")
    if not args.no_sbml:
        print(f"[SUCCESS] Model SBML disimpan ke: {args.sbml}")

    if args.sensitivity:
        payload = run_sensitivity(
            glucose_limit=args.glucose_limit,
            oxygen_limit=args.oxygen_limit,
            tyrosine_limit=args.tyrosine_limit,
            atp_maintenance=args.atp_maintenance,
            biomass_requirement=args.biomass,
        )
        print(f"[SUCCESS] Sensitivity analysis disimpan ke: {SENSITIVITY_RESULTS_PATH}")
        print(f"          {len(payload['rows'])} kombinasi diuji.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
