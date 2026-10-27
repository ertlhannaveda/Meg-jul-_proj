"""Nevleges, szubkritikus szerves Rankine-ciklus (ORC) modell a CoolProp hasznalataval.

A modul hatokore: termodinamikai allapotok, ciklusteljesitmeny, az egyes komponensek izentropikus hatasfoka,
valamint az elparologtato pinch-pontjanak ellenorzese. Ez a modul NEM vegzi el a 8760 oras idoszakra vonatkozo modellezest,
az elparologtatasi nyomas optimalizalasat, az adott evre legmegfelelobb munkakozeg kivalasztasat,
sem pedig a kut-szivattyuk vagy a szarazhuto-ventilatorok energiaigenyenek levonasat
ezek az eves modell feladatkorebe tartoznak.

A mertekegysegek minden visszaadott ertek eseteben fel vannak tuntetve.
"""

from __future__ import annotations

from math import isfinite

import numpy as np
from CoolProp.CoolProp import PropsSI


DEFAULT_FLUIDS = ("R245fa", "R1233zd(E)", "R600a")


def evaporation_pressure_bar(fluid: str, T_evap_C: float) -> float:
    """Az adott munkakozeg parolgasi nyomasanak meghatarozasa a parolgasi homerseklet alapjan.

    Lehetove teszi a kulonbozo munkakozegek osszehasonlitasat
    azonos nevleges parolgasi homersekleten.
    Hanna a kesobbi optimalizalas soran tetszoleges, fizikailag
    megengedett parolgasi nyomast (p_evap_bar) is megadhat
    a calculate_orc fuggvenynek.
    """
    return PropsSI("P", "T", T_evap_C + 273.15, "Q", 0, fluid) / 1e5


def _state(fluid: str, pressure_Pa: float, enthalpy_J_kg: float) -> dict[str, float]:
    """termodinamikai allapotok mernoki mertekegysegekben"""
    return {
        "T_C": PropsSI("T", "P", pressure_Pa, "H", enthalpy_J_kg, fluid) - 273.15,
        "p_bar": pressure_Pa / 1e5,
        "h_kJ_kg": enthalpy_J_kg / 1e3,
        "s_kJ_kgK": PropsSI("S", "P", pressure_Pa, "H", enthalpy_J_kg, fluid) / 1e3,
    }


def calculate_orc(
    *,
    fluid: str,
    T_source_C: float,
    T_air_C: float,
    m_brine_kg_s: float,
    p_evap_bar: float,
    eta_turbine: float = 0.80,
    eta_pump: float = 0.70,
    pinch_K: float = 8.0,
    superheat_K: float = 3.0,
    condenser_approach_K: float = 15.0,
    min_condenser_C: float = 15.0,
    T_brine_out_min_C: float = 65.0,
    cp_brine_J_kgK: float = 4180.0,
) -> dict:
    """Szamitson ki EGY allandosult allapotu ORC uzemi pontot termelo kutankent.

            1: kondenzator kimenet, telitett folyadek p_cond nyomason
            2: ORC betaplalo szivattyu kimenet, tulhutott suritett folyadek p_evap homersekleten
            3: elparologtato kimenet / turbina bemenet, tulhevitett goz p_evap homersekleten
            4: valos turbina kimenet p_cond homersekleten

            A szerves folyadek tomegaramat a minimalis sooldat-kimeneti homerseklet es az ellenaramu HX szukuleti korlatok kozul
            a kisebbkent valasztjuk.

            A sooldat egy kulonallo nyomas alatti vizkor (allando cp modell).

            A p_evap_bar egy BEMENET, itt nem optimalizalva. A visszatero ciklusteljesitmeny a turbina tengelyteljesitmenye
            minusz a belso ORC betaplalo szivattyu teljesitmenye.

            A generator vesztesegeit, a sooldat-kut szivattyuit, a hutoventilatorokat es az eves rendelkezesre allast
            szandekosan a downstream rendszermodellhez hagytuk.
    """
    inputs = [T_source_C, T_air_C, m_brine_kg_s, p_evap_bar, eta_turbine,
              eta_pump, pinch_K, superheat_K, condenser_approach_K,
              min_condenser_C, T_brine_out_min_C, cp_brine_J_kgK]
    if not all(isfinite(float(value)) for value in inputs):
        raise ValueError("All numerical inputs must be finite.")
    if m_brine_kg_s <= 0 or cp_brine_J_kgK <= 0 or p_evap_bar <= 0:
        raise ValueError("Mass flow, heat capacity and evaporation pressure must be positive.")
    if not (0 < eta_turbine <= 1 and 0 < eta_pump <= 1):
        raise ValueError("Isentropic component efficiencies must lie in (0, 1].")
    if min(pinch_K, superheat_K, condenser_approach_K) < 0:
        raise ValueError("Pinch, superheat and condenser approach cannot be negative.")
    if T_brine_out_min_C >= T_source_C:
        raise ValueError("Brine outlet minimum must be lower than source temperature.")

    p_evap_Pa = p_evap_bar * 1e5
    T_cond_C = max(min_condenser_C, T_air_C + condenser_approach_K)
    p_cond_Pa = PropsSI("P", "T", T_cond_C + 273.15, "Q", 0, fluid)
    p_crit_Pa = PropsSI("PCRIT", fluid)
    if not (p_cond_Pa < p_evap_Pa < p_crit_Pa):
        raise ValueError("A subcritical cycle requires p_cond < p_evap < p_crit.")

    T_evap_C = PropsSI("T", "P", p_evap_Pa, "Q", 0, fluid) - 273.15
    T_in_turbine_C = T_evap_C + superheat_K
    if T_source_C - T_in_turbine_C < pinch_K:
        raise ValueError("The hot end of the evaporator violates the pinch constraint.")

    # 1 -> 2: ORC tapszivattyu izentropikus referencia vegponttal 2s
    h1 = PropsSI("H", "P", p_cond_Pa, "Q", 0, fluid)
    s1 = PropsSI("S", "P", p_cond_Pa, "Q", 0, fluid)
    h2s = PropsSI("H", "P", p_evap_Pa, "S", s1, fluid)
    h2 = h1 + (h2s - h1) / eta_pump

    # 2 -> 3: elomelegites allando nyomason, forras es meghatarozott tulhevites
    h3 = PropsSI("H", "P", p_evap_Pa, "T", T_in_turbine_C + 273.15, fluid)
    s3 = PropsSI("S", "P", p_evap_Pa, "H", h3, fluid)

    # 3 -> 4: turbina-tagulas, izentropikus 4s referencia-vegallapottal.
    h4s = PropsSI("H", "P", p_cond_Pa, "S", s3, fluid)
    h4 = h3 - eta_turbine * (h3 - h4s)

    h_sat_liq = PropsSI("H", "P", p_evap_Pa, "Q", 0, fluid)
    h_sat_vap = PropsSI("H", "P", p_evap_Pa, "Q", 1, fluid)
    if not (h1 < h2 < h_sat_liq < h_sat_vap < h3 and h4 > h1):
        raise ValueError("Unexpected state order for the assumed simple subcritical ORC.")

    # Ellenaramu hocserelo: a meleg kozeg a 3. allapothoz kozel lep be.
    # A munkakozeg h entalpiajanal, a nagynyomasu oldalon:
    # T_brine(h) = T_source - m_orc*(h3-h)/(m_brine*c_brine).
    # A T_brine(h)-T_ORC(h) >= pinch feltetel meghatarozza az m_orc felso korlatjat.
    enthalpy_grid = np.unique(np.concatenate([
        np.linspace(h2, h_sat_liq, 16),
        np.linspace(h_sat_liq, h_sat_vap, 41),
        np.linspace(h_sat_vap, h3, 16),
    ]))
    heat_capacity_rate_W_K = m_brine_kg_s * cp_brine_J_kgK
    pinch_flow_limits = []
    for enthalpy in enthalpy_grid:
        h_drop = h3 - enthalpy
        if h_drop < 1e-7:
            continue
        T_fluid_C = PropsSI("T", "P", p_evap_Pa, "H", float(enthalpy), fluid) - 273.15
        numerator = heat_capacity_rate_W_K * (T_source_C - pinch_K - T_fluid_C)
        if numerator <= 0:
            raise ValueError("No positive ORC mass flow can satisfy the HX pinch.")
        pinch_flow_limits.append(numerator / h_drop)

    brine_out_flow_limit = (heat_capacity_rate_W_K * (T_source_C - T_brine_out_min_C)
                            / (h3 - h2))
    m_orc_kg_s = float(min(min(pinch_flow_limits), brine_out_flow_limit))
    if m_orc_kg_s <= 0:
        raise ValueError("Nonpositive calculated working-fluid flow.")

    T_brine_out_C = T_source_C - m_orc_kg_s * (h3 - h2) / heat_capacity_rate_W_K
    hx_temperature_differences = []
    for enthalpy in enthalpy_grid:
        T_fluid_C = PropsSI("T", "P", p_evap_Pa, "H", float(enthalpy), fluid) - 273.15
        T_brine_C = T_source_C - m_orc_kg_s * (h3 - enthalpy) / heat_capacity_rate_W_K
        hx_temperature_differences.append(T_brine_C - T_fluid_C)
    min_actual_pinch_K = float(min(hx_temperature_differences))
    if min_actual_pinch_K < pinch_K - 1e-4:
        raise AssertionError("Numerical pinch constraint not satisfied.")

    # Az elso fotetel szerinti egyensuly: q_evap - q_cond = w_turbine - w_feed_pump.
    Q_evap_W = m_orc_kg_s * (h3 - h2)
    P_turbine_W = m_orc_kg_s * (h3 - h4)
    P_orc_pump_W = m_orc_kg_s * (h2 - h1)
    Q_cond_W = m_orc_kg_s * (h4 - h1)
    P_cycle_net_W = P_turbine_W - P_orc_pump_W
    if abs(Q_evap_W - Q_cond_W - P_cycle_net_W) > 1e-7 * Q_evap_W:
        raise AssertionError("The first-law energy balance does not close.")
    if P_cycle_net_W <= 0 or Q_cond_W <= 0:
        raise ValueError("The operating point has no usable positive cycle power.")

    # A goztartalomnak csak a ketfazisu tartomanyon belul van ertelme.
    h_cond_sat_vap = PropsSI("H", "P", p_cond_Pa, "Q", 1, fluid)
    x4 = None
    if h1 < h4 < h_cond_sat_vap:
        x4 = float((h4 - h1) / (h_cond_sat_vap - h1))

    states = {
        # Telitett allapot: a P,Q-bol szarmaztatott ertakeket kozvetlenul kell megadni,
        # elkerulve ezzel a P,H-villanas (flash) ketertelmuseget pontosan a telitesi hatarvonalon.
        "1": {
            "T_C": T_cond_C,
            "p_bar": p_cond_Pa / 1e5,
            "h_kJ_kg": h1 / 1e3,
            "s_kJ_kgK": s1 / 1e3,
        },
        "2": _state(fluid, p_evap_Pa, h2),
        "3": _state(fluid, p_evap_Pa, h3),
        "4": _state(fluid, p_cond_Pa, h4),
    }
    return {
        "fluid": fluid,
        "states": states,
        "isentropic_reference": {
            "2s": _state(fluid, p_evap_Pa, h2s),
            "4s": _state(fluid, p_cond_Pa, h4s),
        },
        "performance": {
            "m_orc_kg_s": m_orc_kg_s,
            "Q_evaporator_kW": Q_evap_W / 1e3,
            "Q_condenser_kW": Q_cond_W / 1e3,
            "P_turbine_kW": P_turbine_W / 1e3,
            "P_orc_pump_kW": P_orc_pump_W / 1e3,
            "P_cycle_net_kW": P_cycle_net_W / 1e3,
            "eta_cycle": P_cycle_net_W / Q_evap_W,
            "T_evap_C": T_evap_C,
            "T_cond_C": T_cond_C,
            "T_brine_out_C": T_brine_out_C,
            "min_actual_pinch_K": min_actual_pinch_K,
            "turbine_exit_quality": x4,
            "energy_balance_residual_kW": (Q_evap_W - Q_cond_W - P_cycle_net_W) / 1e3,
        },
        "assumptions": {
            "T_source_C": T_source_C,
            "T_air_C": T_air_C,
            "m_brine_kg_s": m_brine_kg_s,
            "p_evap_bar": p_evap_bar,
            "eta_turbine": eta_turbine,
            "eta_pump": eta_pump,
            "pinch_K": pinch_K,
            "superheat_K": superheat_K,
            "condenser_approach_K": condenser_approach_K,
            "T_brine_out_min_C": T_brine_out_min_C,
            "cp_brine_J_kgK": cp_brine_J_kgK,
            "auxiliaries_included": "ORC feed pump only; no condenser fan or geothermal well pumps",
        },
    }


if __name__ == "__main__":
    print("NOMINAL FLUID SCREENING (fixed saturation temperature 85 C)")
    for candidate in DEFAULT_FLUIDS:
        pressure = evaporation_pressure_bar(candidate, 85.0)
        try:
            case = calculate_orc(fluid=candidate, T_source_C=115.0,
                                 T_air_C=13.53, m_brine_kg_s=40.0,
                                 p_evap_bar=pressure)
            perf = case["performance"]
            print(f"{candidate:12s} p_evap={pressure:7.3f} bar  "
                  f"P_cycle={perf['P_cycle_net_kW']:8.1f} kW  "
                  f"eta={100*perf['eta_cycle']:5.2f}%  "
                  f"pinch={perf['min_actual_pinch_K']:5.2f} K")
        except ValueError as exc:
            print(f"{candidate:12s} NOT FEASIBLE: {exc}")
