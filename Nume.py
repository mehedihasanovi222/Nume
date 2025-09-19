import math

# Step size
h = 0.1
# Number of steps
N = 5

x_val = [0]
y_euler = [1]
y_rk4 = [1]
y_exact = [1]
abs_error_euler = [0]
abs_error_rk4 = [0]
rel_error_euler = [0]
rel_error_rk4 = [0]
accuracy_note_list = []
stability_note_list = []

# ODE function
def f(x, y):
    """
    ODE function for dy/dx = -2y
    Parameters:
        x (float)
        y (float)
    Returns:
        float
    """
    return -2 * y

# Compute step by step
for i in range(N):
    xn = x_val[-1]
    yn_e = y_euler[-1]
    yn_r = y_rk4[-1]

    # Euler method
    y_next_e = yn_e + h * f(xn, yn_e)
    y_euler.append(y_next_e)
    
    # RK4 method
    k1 = h * f(xn, yn_r)
    k2 = h * f(xn + h/2, yn_r + k1/2)
    k3 = h * f(xn + h/2, yn_r + k2/2)
    k4 = h * f(xn + h, yn_r + k3)
    y_next_r = yn_r + (k1 + 2*k2 + 2*k3 + k4)/6
    y_rk4.append(y_next_r)

    # Exact solution
    x_next = xn + h
    y_ex = math.exp(-2 * x_next)
    y_exact.append(y_ex)
    
    # Error Euler and RK4
    abs_error_euler.append(abs(y_next_e - y_ex))
    abs_error_rk4.append(abs(y_next_r - y_ex))
    rel_error_euler.append(abs(y_next_e - y_ex) / y_ex)
    rel_error_rk4.append(abs(y_next_r - y_ex) / y_ex)    

    # Accuracy and stability notes for this step
    accuracy_note_list.append('Euler: less accurate; RK4: highly accurate')
    stability_note_list.append('Euler stable if h < 1; RK4 more stable')
    
    # Update x
    x_val.append(x_next)

# Print Table Header
print(f"{'x':<6}{'Exact y':<12}{'Euler y':<12}{'Euler Abs Err':<15}{'Euler Rel Err':<15}"
      f"{'RK4 y':<12}{'RK4 Abs Err':<15}{'RK4 Rel Err':<15}{'Notes':<40}")

# Print Table Row
for i in range(N+1):
    print(f"{x_val[i]:<6.2f}{y_exact[i]:<12.5f}{y_euler[i]:<12.5f}"
          f"{abs_error_euler[i]:<15.5f}{rel_error_euler[i]:<15.5f}"
          f"{y_rk4[i]:<12.5f}{abs_error_rk4[i]:<15.5f}{rel_error_rk4[i]:<15.5f}"
          f"{accuracy_note_list[i-1] if i > 0 else 'Euler: N/A; RK4: N/A'} | {stability_note_list[i-1] if i > 0 else 'Euler: N/A; RK4: N/A'}")

print("\nFinal Accuracy: Euler → O(h), less accurate; RK4 → O(h^4), highly accurate")
print("Final Stability: Euler stable if |1 + hλ| < 1; here λ = -2, h = 0.1 → Stable; RK4 more stable overall")
