"""Tracked complete-quartic reconstruction and finite-mode one-loop two-point.

The physical input is the same ten-mode reduced Gaunt Hamiltonian. The result
is a bare, Weyl-prescription one-loop two-point coefficient, not a four-graviton
amplitude or a proof of interacting reconstruction.
"""

from fractions import Fraction as F
from math import pi, sqrt
from pathlib import Path
import json
import os
import platform
import sys
import time

from run import report
from complete_quartic import (
    calculate_coefficients, complete_polynomials, canonical_factors,
    hessian_block, json_coefficients, matrix_key, lt, lm,
)
from reduced_gaunt_exact import five_polarizations

ROOT=Path(__file__).resolve().parent


def scalar_number(value):
    value=complex(value)
    if abs(value.imag)>1e-10*(1+abs(value.real)):
        raise ArithmeticError("Unexpected imaginary Euclidean two-point coefficient")
    return value.real


def analytic_two_point(frequencies, block1, block2, quartic, omega):
    tadpole=sum(quartic[2*i][2*i]/(4*w)+w*quartic[2*i+1][2*i+1]/4
                for i,w in enumerate(frequencies))
    bubble=0.0
    for i,a in enumerate(frequencies):
        for j,b in enumerate(frequencies):
            integral=(a+b)/(2*a*b*((a+b)**2+omega**2))
            bubble-=block1[2*i][2*j]*block2[2*j][2*i]*integral/2
    # All gravitational cubic tensors here are pure configuration tensors.
    assert all(not block1[i][j] and not block2[i][j]
               for i in range(20) for j in range(20) if i%2 or j%2)
    return {"tadpole":tadpole,"bubble":bubble,"total":tadpole+bubble}


def polynomial_historical_check(cubic,quartic):
    old=json.loads((ROOT.parents[1]/"reports/quantum-bridge/data/reduced-vertices-results.json").read_text())
    basis=five_polarizations()
    output=[]
    for row in old["exact_directions"]:
        labels=[matrix_key(a) for a in row["labels_Eq_Bq_Es_Bs"]]
        variables=[F(0)]*20
        for parity in (0,1):
            for i,b in enumerate(basis):
                norm=lt(lm(b,b))
                variables[2*(5*parity+i)]=lt(lm(labels[parity],b))/norm
                variables[2*(5*parity+i)+1]=lt(lm(labels[2+parity],b))/norm
        def evaluate(poly):
            value=F(0)
            for key,coefficient in poly.items():
                for index in key:
                    coefficient*=variables[index]
                value+=coefficient
            return value
        assert evaluate(cubic)==F(row["coefficients"]["H3_configuration"])
        target=sum((F(row["coefficients"][key]) for key in
                    ("H4_configuration","H4_mixed","H4_momentum")),F(0))
        assert evaluate(quartic)==target,row["name"]
        output.append({"direction":row["name"],"full_polynomial_residual":"0"})
    return output


def main(config):
    from phase_space_loop import integrate_vertex
    from loop_kernel_controls import run_controls
    start=time.monotonic()
    print(json.dumps({"configuration":config,"python":sys.version,"platform":platform.platform(),
        "commit":os.getenv("ORX_COMMIT"),"run":os.getenv("ORX_RUN_ID"),
        "experiment":os.getenv("ORX_EXPERIMENT")}),flush=True)
    report(phase="quartic_coefficients",step=0,of=24,unit="directions")
    counter=0
    def progress(row):
        nonlocal counter
        counter+=1
        print("EXACT QUARTIC DIRECTION "+json.dumps(row),flush=True)
        report(step=counter,m={"coefficient_directions":counter})
    coefficients,raw,old_checks=calculate_coefficients(progress)
    assert counter == 24
    print("COMPLETE INVARIANT COEFFICIENTS "+json.dumps(json_coefficients(coefficients)),flush=True)
    report(phase="canonical_polynomial",step=0,of=1,unit="checks")
    cubic,quartic=complete_polynomials(coefficients)
    polynomial_checks=polynomial_historical_check(cubic,quartic)
    gram,factors,frequencies=canonical_factors()
    unit3=(2*pi/3)**1.5
    unit4=pi**1.5/sqrt(2)
    one_point=[]
    for external in range(20):
        b=hessian_block(cubic,(external,),unit3,factors)
        value=sum(b[2*i][2*i]/(4*w)+w*b[2*i+1][2*i+1]/4 for i,w in enumerate(frequencies))
        one_point.append(value)
    assert max(abs(x) for x in one_point)<config["loop_relative_tolerance"]
    report(step=1,m={"quartic_monomials":len(quartic),"cubic_monomials":len(cubic),
                      "exact_polynomial_residual":0})
    print("CANONICAL POLYNOMIAL "+json.dumps({"cubic_monomials":len(cubic),
        "quartic_monomials":len(quartic),"historical_checks":polynomial_checks,
        "maximum_one_point_coefficient":max(abs(x) for x in one_point)}),flush=True)

    report(phase="loop_kernel_controls",step=0,of=1,unit="suites")
    controls=run_controls(config)
    print("LOOP KERNEL CONTROLS "+json.dumps(controls),flush=True)
    if not controls["passed"]:
        raise ArithmeticError("One-loop bookkeeping controls failed; see the full control record above")
    report(step=1,m={"loop_control_suites":1})

    report(phase="gravitational_two_point",step=0,of=22,unit="cases")
    # Six frequency points per parity (Q and P) plus four rotational
    # comparisons per parity. All free internal modes participate.
    specs=[]
    for mode in (1,6):
        for omega in config["external_frequencies"]:
            for kind in (0,1):
                specs.append((mode,kind,omega))
    for mode in (0,3,5,8):
        for kind in (0,1):
            specs.append((mode,kind,config["external_frequencies"][1]))
    assert len(specs)==20
    rows=[]
    for count,(mode,kind,omega) in enumerate(specs,1):
        index=2*mode+kind
        b=hessian_block(cubic,(index,),unit3,factors)
        c=hessian_block(quartic,(index,index),unit4,factors)
        blocks={frozenset((0,)):b,frozenset((1,)):b,frozenset((0,1)):c}
        target=analytic_two_point(frequencies,b,b,c,omega)
        estimates=[scalar_number(integrate_vertex(frequencies,[omega,-omega],blocks,quadrature_order=n))
                   for n in config["loop_quadrature_orders"]]
        error=abs(estimates[-1]-target["total"])/(1+abs(target["total"]))
        assert error<config["loop_relative_tolerance"],(mode,kind,omega,error)
        row={"mode":mode,"external_variable":"Q" if kind==0 else "P","external_frequency":omega,
             "tensor_parity":"electric" if mode<5 else "magnetic",**target,
             "quadrature_orders":config["loop_quadrature_orders"],"quadrature_values":estimates,
             "quadrature_difference":abs(estimates[-1]-estimates[0]),"scaled_error":error}
        rows.append(row)
        print("GRAVITATIONAL TWO POINT "+json.dumps(row),flush=True)
        report(row="two_point",m={"case":f"{mode}:{row['external_variable']}:{omega}",
            "external_frequency":omega,"self_energy":target["total"],"loop_error":error})
        report(step=count,m={"loop_error":error,"gravitational_cases":count})
    # Rotation-equivariance predicts equal normalized results within each
    # irreducible multiplet, including its different magnetic components.
    rotational_errors=[]
    for row in rows:
        if row["mode"] in (1,6):
            continue
        representative=1 if row["mode"]<5 else 6
        reference=next(x for x in rows if x["mode"]==representative and
                       x["external_variable"]==row["external_variable"] and
                       x["external_frequency"]==row["external_frequency"])
        error=abs(row["total"]-reference["total"])/(1+abs(reference["total"]))
        rotational_errors.append(error)
        assert error<config["loop_relative_tolerance"],(row,reference)
    mixed_phase=[]
    for count,mode in enumerate((1,6),21):
        omega=config["external_frequencies"][1]
        bq=hessian_block(cubic,(2*mode,),unit3,factors)
        bp=hessian_block(cubic,(2*mode+1,),unit3,factors)
        bqp=hessian_block(quartic,(2*mode,2*mode+1),unit4,factors)
        blocks={frozenset((0,)):bq,frozenset((1,)):bp,frozenset((0,1)):bqp}
        value=scalar_number(integrate_vertex(frequencies,[omega,-omega],blocks,
                                           quadrature_order=config["loop_quadrature_orders"][-1]))
        assert abs(value)<config["loop_relative_tolerance"]
        row={"mode":mode,"external_frequency":omega,"Pi_QP":value,"predicted_Pi_QP":0}
        mixed_phase.append(row)
        print("MIXED PHASE CHECK "+json.dumps(row),flush=True)
        report(step=count,m={"mixed_phase_error":abs(value)})
    configuration=[]
    for mode in (1,6):
        for omega in config["external_frequencies"]:
            selected={row["external_variable"]:row for row in rows
                      if row["mode"]==mode and row["external_frequency"]==omega}
            configuration.append({"mode":mode,"external_frequency":omega,
                "Q_inverse_propagator_correction":selected["Q"]["total"]-omega**2*selected["P"]["total"],
                "formula":"Pi_QQ(Omega) - Omega^2 Pi_PP; hbar coefficient in normalized canonical variables"})

    report(phase="export",step=0,of=1,unit="reports")
    maximum=max(row["scaled_error"] for row in rows)
    result={"config":config,"source_commit":os.getenv("ORX_COMMIT"),"run":os.getenv("ORX_RUN_ID"),
        "physical_free_rank_N":3,"selected_free_modes":10,"free_frequencies":frequencies,
        "gram_diagonal":gram,"canonical_field_factors":factors,
        "invariant_coefficients":json_coefficients(coefficients),"exact_directions":raw,
        "historical_functional_checks":old_checks,"historical_polynomial_checks":polynomial_checks,
        "cubic_raw_polynomial":[{"indices":list(k),"coefficient":str(v)} for k,v in sorted(cubic.items())],
        "quartic_raw_polynomial":[{"indices":list(k),"coefficient":str(v)} for k,v in sorted(quartic.items())],
        "polynomial_units":{"cubic":"(2*pi/3)^(3/2)","quartic":"pi^(3/2)/sqrt(2)",
                            "variables":"interleaved unnormalized q_alpha and S_alpha; multiply the listed canonical factors"},
        "loop_kernel_controls":controls,"gravitational_two_point":rows,
        "mixed_phase_checks":mixed_phase,"configuration_inverse_propagator":configuration,
        "maximum_loop_scaled_error":maximum,"maximum_rotational_scaled_error":max(rotational_errors),
        "quantization":"Bare midpoint/Weyl perturbative prescription; hbar coefficient; no order-hbar counterterm added",
        "limitations":[
            "Finite ten-mode two-point function, not a four-graviton scattering amplitude.",
            "All four-point topologies are tested on controls, not yet evaluated with complete gravitational H5/H6.",
            "One radial profile per parity; no interacting physical large-N limit is established.",
            "Euclidean frequency-dependent coefficients do not prove a positive full interacting Hamiltonian or reconstruction.",
            "Finite-basis frequency shifts must not be identified with a graviton mass.",
            "Full local corner quantization and microscopic state counting remain open."
        ],"status":"all specified complete-quartic and finite-mode loop checks passed",
        "wall_seconds":time.monotonic()-start}
    output=ROOT/"output"
    output.mkdir(exist_ok=True)
    path=output/"quartic_loop_results.json"
    path.write_text(json.dumps(result,indent=2)+"\n")
    report(step=1,m={"maximum_loop_error":maximum,"maximum_rotation_error":max(rotational_errors),
                      "quartic_invariant_coefficients":27})
    print("FINAL "+json.dumps({"result":str(path),"source_commit":result["source_commit"],
        "quartic_invariant_coefficients":27,"cubic_monomials":len(cubic),"quartic_monomials":len(quartic),
        "gravitational_two_point_cases":len(rows),"maximum_loop_scaled_error":maximum,
        "maximum_rotational_scaled_error":max(rotational_errors),"status":result["status"],
        "scope":"Complete finite-mode one-loop two-point; four-graviton amplitude remains uncomputed",
        "wall_seconds":result["wall_seconds"]}),flush=True)
