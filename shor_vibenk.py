# pylint: skip-file

#Shor-Kitaev próbaalkalmazás


#A külső, nyilvánosan elérhető könyvtárak betöltése
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


# A kvantumalgoritmus
def shor_kvantum(eng, N, a, verbose=False):
    """
    Az algoritmus kvantum része

    Param:
        eng (MainEngine): A fő engine
        N (int): A faktorizálandó egész szán
        a (int): A próbálkozás bázisaként használt relatív prím
        verbose (bool): A köztes mérési eredmények mutatása?

    Visszatér:
        r (float): Az a bázison alapuló periodicitás
    """
    n = int(math.ceil(math.log(N, 2)))
    print("A szükséges kvantumregiszter mérete: ",n," kvantumbit")

    x = eng.allocate_qureg(n)

    X | x[0]       # Pauli-X kapu

    measurements = [0] * (2 * n) # a mérési eredmények tárolására

    ctrl_qubit = eng.allocate_qubit()

    for k in range(2 * n):  # egy QPE iteráció
        current_a = pow(a, 1 << (2 * n - 1 - k), N)
        H | ctrl_qubit      # Hadamard kapu
        with Control(eng, ctrl_qubit):
            MultiplyByConstantModN(current_a, N) | x

        # QFT --> a mérési eredmények függvényében forgatás
        for i in range(k):
            if measurements[i]:
                R(-math.pi / (1 << (k - i))) | ctrl_qubit
        H | ctrl_qubit      # Hadamard kapu

        Measure | ctrl_qubit
        eng.flush()
        measurements[k] = int(ctrl_qubit)
        if measurements[k]:
            X | ctrl_qubit      # Pauli-X kapu

        if verbose:
            print(f"\033[95m{measurements[k]}\033[0m", end="")
            sys.stdout.flush()

    All(Measure) | x    # mérési eredmények [0,1)
    y = sum((measurements[2 * n - 1 - i] * 1.0 / (1 << (i + 1))) for i in range(2 * n))

    # LNKO
    r = Fraction(y).limit_denominator(N - 1).denominator

    # a lehetséges periódus kiírása
    return r


def high_level_gates(eng, cmd):
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

    # a szimuláció létrehozása és paraméterezése
    eng = MainEngine(Simulator(), compilerengines)

    # futási paraméterek bekérése
    N = int(input('\n\tA faktorializálandó szám (N): '))

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
        # kvantum rész futtatása
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
        else:
            print(f"\n\n\t\033[91mAz algoritmus nem jutott eredményre a megadott próbák számán belül! Kérjük próbálja másik kiinduló tippel!")