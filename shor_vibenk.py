# pylint: skip-file

"""Example implementation of Shor's algorithm."""

import math
import random
import sys
from fractions import Fraction

try:
    from math import gcd
except ImportError:
    from fractions import gcd

import projectq.libs.math
import projectq.setups.decompositions
from projectq.backends import ResourceCounter, Simulator
from projectq.cengines import (
    AutoReplacer,
    DecompositionRuleSet,
    InstructionFilter,
    LocalOptimizer,
    MainEngine,
    TagRemover,
)
from projectq.libs.math import AddConstant, AddConstantModN, MultiplyByConstantModN
from projectq.meta import Control
from projectq.ops import QFT, All, BasicMathGate, H, Measure, R, Swap, X, get_inverse


def shor_kvantum(eng, N, a, verbose=False):
    """
    Run the quantum subroutine of Shor's algorithm for factoring.

    Args:
        eng (MainEngine): Main compiler engine to use.
        N (int): Number to factor.
        a (int): Relative prime to use as a base for a^x mod N.
        verbose (bool): If True, display intermediate measurement results.

    Returns:
        r (float): Potential period of a.
    """
    n = int(math.ceil(math.log(N, 2)))
    print("A szükséges kvantumregiszter mérete: ",n," kvantumbit")

    x = eng.allocate_qureg(n)

    X | x[0]

    measurements = [0] * (2 * n)  # will hold the 2n measurement results

    ctrl_qubit = eng.allocate_qubit()

    for k in range(2 * n):
        current_a = pow(a, 1 << (2 * n - 1 - k), N)
        # one iteration of 1-qubit QPE
        H | ctrl_qubit
        with Control(eng, ctrl_qubit):
            MultiplyByConstantModN(current_a, N) | x

        # perform inverse QFT --> Rotations conditioned on previous outcomes
        for i in range(k):
            if measurements[i]:
                R(-math.pi / (1 << (k - i))) | ctrl_qubit
        H | ctrl_qubit

        # and measure
        Measure | ctrl_qubit
        eng.flush()
        measurements[k] = int(ctrl_qubit)
        if measurements[k]:
            X | ctrl_qubit

        if verbose:
            print(f"\033[95m{measurements[k]}\033[0m", end="")
            sys.stdout.flush()

    All(Measure) | x
    # turn the measured values into a number in [0,1)
    y = sum((measurements[2 * n - 1 - i] * 1.0 / (1 << (i + 1))) for i in range(2 * n))

    # continued fraction expansion to get denominator (the period?)
    r = Fraction(y).limit_denominator(N - 1).denominator

    # return the (potential) period
    return r


# Filter function, which defines the gate set for the first optimization
# (don't decompose QFTs and iQFTs to make cancellation easier)
def high_level_gates(eng, cmd):
    """Filter high-level gates."""
    g = cmd.gate
    if g == QFT or get_inverse(g) == QFT or g == Swap:
        return True
    if isinstance(g, BasicMathGate):
        return False
        if isinstance(g, AddConstant):
            return True
        elif isinstance(g, AddConstantModN):
            return True
        return False
    return eng.next_engine.is_available(cmd)


if __name__ == "__main__":
    # build compilation engine list
    resource_counter = ResourceCounter()
    rule_set = DecompositionRuleSet(modules=[projectq.libs.math, projectq.setups.decompositions])
    compilerengines = [
        AutoReplacer(rule_set),
        InstructionFilter(high_level_gates),
        TagRemover(),
        LocalOptimizer(3),
        AutoReplacer(rule_set),
        TagRemover(),
        LocalOptimizer(3),
        resource_counter,
    ]

    # make the compiler and run the circuit on the simulator backend
    eng = MainEngine(Simulator(), compilerengines)

    # print welcome message and ask the user for the number to factor
    N = int(input('\n\tA faktorializálandó szám (N): '))
#    print(f"\n\tFactoring N = {N}: \033[0m", end="")

    # choose a base at random:
    g = int(input("A választott kiinduló tipp (g) 0=RNG: "))
    if g==0:
        g = int(random.random() * N)

    print("\n\n***** Faktorializálási paraméterek *****")
    print("  N =",N)
    print("  g =",g)
    print("  r = ?")


    if not gcd(g, N) == 1:
        print("\n\n\t\033[92mNagy szerencse, a kiinduló tipp már az egyik faktor :)")
    else:
        # run the quantum subroutine

        kor = 0
        r = -1
        print("\n\n***** Kvantum algoritmus indítása *****")

        while r % 2 != 0 or r < 3 and kor < 6:
            kor = kor + 1
            print("\n\n\tPróbálkozás száma : #",kor)
            r = shor_kvantum(eng, N, g, True)
            print("\nMért periódus: r=",r)

        print("\n\n***** Futási eredmény kiértékelése *****")
        print("\nPeriodicitás (r) értéke=", r)


        
        # try to determine the factors
#        if r % 2 == 0 and r>2:
        if kor < 6:
            tippalap = pow(g, r >> 1)

            tipp1 = tippalap + 1
            tipp2 = tippalap - 1

            print("A két jobb tipp: f1=",tipp1,"; f2=",tipp2)
            f1 = gcd(tipp1, N)
            print(tipp1," és ", N, "legnagyobb közös osztója : ",f1)
            f2 = gcd(tipp2, N)
            print(tipp2," és ", N, "legnagyobb közös osztója : ",f2)

            if (not f1 * f2 == N) and f1 * f2 > 1 and int(1.0 * N / (f1 * f2)) * f1 * f2 == N:
                f1, f2 = f1 * f2, int(N / (f1 * f2))
            if f1 * f2 == N and f1 > 1 and f2 > 1:
                print(f"\n\n\t\033[92mPrímtényezők megtalálva : {f1} * {f2} = {N}\033[0m")
            else:
                print(f"\n\n\t\033[91mKérjük próbálja másik kiinduló tippel!")

            #print(resource_counter)  # print resource usage
        else:
            print(f"\n\n\t\033[91mAz algoritmus nem jutott eredményre a megadott próbák számán belül! Kérjük próbálja másik kiinduló tippel!")