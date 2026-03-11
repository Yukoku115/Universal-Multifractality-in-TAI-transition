import numpy as np

data      = np.load("phase_diagram_checkpoint.npz")
done_mask = data["done_mask"]
M_vals    = data["M_vals"]
W_vals    = data["W_vals"]

completed = int(done_mask.sum())
total     = done_mask.size
print(f"Completed: {completed}/{total} ({100*completed/total:.1f}%)")

for i, m in enumerate(M_vals):
    n = done_mask[i].sum()
    if n == 0:
        status = "not started"
    elif n == len(W_vals):
        status = "complete"
    else:
        status = f"partial — {n}/{len(W_vals)} done"
    print(f"  M={m:.3f}  {status}")
