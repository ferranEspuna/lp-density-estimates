import sys
from search_witness_dil import evaluate_bound

def test_aps():
    best_bound = 0.5
    for nx in range(2, 6):
        for ny in range(2, 6):
            if nx + ny > 11: continue
            for dx in range(1, 4):
                for dy in range(1, 4):
                    for start_y in range(0, 4):
                        X = [i * dx for i in range(nx)]
                        Y = [start_y + i * dy for i in range(ny)]
                        
                        has_overlap = False
                        for i in range(nx):
                            for j in range(nx):
                                for k in range(ny):
                                    if X[i] + X[j] == 2 * Y[k]:
                                        has_overlap = True
                                        break
                        if not has_overlap:
                            continue

                        bound = evaluate_bound(X, Y)
                        if bound < best_bound - 1e-4:
                            best_bound = bound
                            print(f"NEW BEST! {best_bound:.5f} with X={X}, Y={Y}", flush=True)

if __name__ == "__main__":
    test_aps()
    print("AP search complete.")
