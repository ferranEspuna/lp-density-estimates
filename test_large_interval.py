import sys
from search_witness_dil import evaluate_bound

def run():
    print("Evaluating large symmetric intervals...")
    for n in range(5, 9):
        X = list(range(n))
        Y = list(range(n))
        print(f"Testing X=[0..{n-1}], Y=[0..{n-1}]")
        bound = evaluate_bound(X, Y)
        print(f"Bound for n={n}: {bound:.5f}", flush=True)

    print("Evaluating disjoint intervals...")
    for n in range(5, 8):
        X = list(range(0, 2*n, 2))
        Y = list(range(1, 2*n, 2))
        print(f"Testing disjoint X={X}, Y={Y}")
        bound = evaluate_bound(X, Y)
        print(f"Bound for disjoint n={n}: {bound:.5f}", flush=True)

if __name__ == "__main__":
    run()
