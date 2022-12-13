"""
https://github.com/plasma-umass/scalene

To profile:
1) pip install scalene
2) cd tests/performance_analysis
3) scalene --profile-all profiling.py
4) cd ...


cProfile:
1) cd tests/performance_analysis
2) python -m cProfile -o profile.stats profiling.py
3) python -m pstats profile.stats
4) cd ..

terminal commands:
1) help
2) sort time   <-- no print
3) stats 10    <-- prints result
https://www.stefaanlippens.net/python_profiling_with_pstats_interactive_mode/

"""

import bigsmiles

polymers = [
    "CC{[>][<]CC(C)[>][<]}CC(C)=C",
    "CCCCCCCCCCCCCCCCCCCCCCCCCCCC",
    "{[>][$]CC[$],[$]CC(CC)[$][<]}",
    "{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}",
    "{[>][<]C(=O)CCCCC(=O)NCCCCCCN[>][<]}",
    "C{[$][$]CC[$],[$]CC(CC)[$][$]}",
    "[H]O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]"
]


def main():

    for i in range(100_000):
        bigsmiles.BigSMILES(polymers[i % len(polymers)])


if __name__ == "__main__":
    main()
