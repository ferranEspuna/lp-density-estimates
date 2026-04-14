import random
from search_witness_dil import evaluate_bound

def main():
    best_bound = 0.5
    size_X = 5
    size_Y = 5
    max_val = 15
    print("Starting randomized search for witness sets...")
    
    attempts = 0
    while True:
        rem_X = random.sample(range(1, max_val + 1), size_X - 1)
        X = [0] + sorted(rem_X)
        Y = sorted(random.sample(range(0, max_val + 1), size_Y))
        
        has_overlap = False
        for i in range(size_X):
            for j in range(size_X):
                for k in range(size_Y):
                    if X[i] + X[j] == 2 * Y[k]:
                        has_overlap = True
                        break
        
        if not has_overlap:
            continue
            
        bound = evaluate_bound(X, Y)
        attempts += 1
        
        if bound < best_bound - 1e-5:
            best_bound = bound
            print(f"[{attempts}] New best bound: {best_bound:.5f} with X={X}, Y={Y}", flush=True)
            if best_bound < 0.49:
                pass
        
        if attempts % 50 == 0:
            print(f"Attempt {attempts}, best so far: {best_bound:.5f}", flush=True)

if __name__ == "__main__":
    main()
