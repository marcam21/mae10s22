import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import brentq

# Inputs
F_ideal = 12.0 # Newtons
P_c_low = 5e5 # Pa
P_c_high = 20e5 # Pa
Pc = np.linspace(P_c_low, P_c_high, 1000)
t_r = 0.005 # seconds

# Constants
T_c = 1222.0 # K
MW = 0.022560 # kg/mol
R = 8.314/MW # J/kg
g = 9.806 # m/s^2
gamma = 1.25
term = (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
C_star = np.sqrt(R*T_c) / (np.sqrt(gamma) * term)

# Given: A_ratios = Ae/At array
A_ratios = np.linspace(30, 100, 1000)

# --- Solve for Mach on both branches (robust bracketing) ---
M_sup = np.zeros_like(A_ratios)
M_sub = np.zeros_like(A_ratios)

for i in range(len(A_ratios)):
    A = A_ratios[i]

    def f(M):
        return (1.0/M) * ( (2.0/(gamma+1)) * (1 + (gamma-1)/2 * M**2) ) ** ( (gamma+1)/(2*(gamma-1)) ) - A

    try:
        M_sub[i] = brentq(f, 0.000001, 1.0 - 1e-9)
    except:
        M_sub[i] = np.nan

    try:
        M_sup[i] = brentq(f, 1.0 + 1e-9, 10000)
    except:
        M_sup[i] = np.nan

# --- Pressure ratios from Mach (isentropic) ---
p_ratio_sup = (1 + (gamma-1)/2 * M_sup**2) ** (-gamma/(gamma-1))
p_ratio_sub = (1 + (gamma-1)/2 * M_sub**2) ** (-gamma/(gamma-1))

CF_ideal = np.sqrt( (2*gamma**2/(gamma-1)) * (2/(gamma+1))**((gamma+1)/(gamma-1)) \
    * (1 - (p_ratio_sup)**((gamma-1)/gamma)) ) \
     + (p_ratio_sup) * A_ratios

u_e = np.sqrt((2*gamma/(gamma-1)*R*T_c)*(1-(p_ratio_sup)**((gamma-1)/gamma)))

T_e = T_c*(p_ratio_sup)**((gamma-1)/gamma)

A_traget = 60.0

CF_ideal_target = interp1d(A_ratios, CF_ideal)(A_traget)
u_e_target = interp1d(A_ratios, u_e)(A_traget)
T_e_target = interp1d(A_ratios, T_e)(A_traget)
p_ratio_target = interp1d(A_ratios, p_ratio_sup)(A_traget)
M_target = interp1d(A_ratios, M_sup)(A_traget)

I_sp_ideal = C_star*CF_ideal_target/g

m_dot_ideal = F_ideal/(I_sp_ideal*g)

A_t = (m_dot_ideal*C_star)/Pc

A_t_target = 4e-6
Pc_target = interp1d(A_t, Pc)(A_t_target)
Pe = (2.0/3.0)*Pc_target*p_ratio_target

V_c = m_dot_ideal * t_r * R * T_c / Pc_target
L_star = V_c/A_t_target
L_c = 0.08
epsilon = L_star/L_c
A_c = A_t_target*epsilon

# sanity check
F = m_dot_ideal*u_e_target + Pc_target*p_ratio_target*A_t_target*A_traget

# Create table
Results = pd.DataFrame({
    'Expansion_Ratio': [A_traget],
    'Chamber_Pressure_Pa': [Pc_target],
    'Exit_Pressure_Pa': [Pe],
    'Exit_Temperature_K': [T_e_target]
})

print(Results)
print("-" * 30)

# --- MoC ---
rt = np.sqrt(A_t_target/np.pi)
N = 48

def moc_minlen_net(gam, Me, Rt, N):

    def nuPM(M):
        val = np.sqrt((gam+1)/(gam-1))*np.arctan( np.sqrt((gam-1)/(gam+1)*(M**2-1)) ) - np.arctan( np.sqrt(M**2-1) )
        return val

    nu_max = 0.5*np.pi*( np.sqrt((gam+1)/(gam-1)) - 1 )

    def invPM(nu_target):
        if not np.isfinite(nu_target) or nu_target <= 1e-9:
            # Handle potential negative nu due to numerical error near 0
            return 1.0
        if nu_target > nu_max:
             # Cap it
             nu_target = nu_max - 1e-9

        def f(M):
            return nuPM(M) - nu_target

        try:
            return brentq(f, 1.0, 1e5)
        except ValueError:
            # If f(b) is still negative, it means nu_target is very close to nu_max
            # and M needs to be huge. Return upper bound.
            if f(1e5) < 0: return 1e5
            # If f(a) is positive (unlikely given checks), return 1.0
            return 1.0

    def muMach(M):
        if M <= 1.0: return np.pi/2.0
        val = 1.0/M
        if val > 1.0: val = 1.0
        return np.arcsin(val)

    # Target turning
    nu_e = nuPM(Me)
    theta_max = 0.5*nu_e
    dtheta = theta_max / N

    # Allocate points structure: List of Lists or 2D Array of Objects/Dicts
    # We need (i, j) where i is Fan Index (1..N), j is Refl Index (1..i)
    # We will use a dictionary mapped by (i, j) tuple
    points = {}

    # Initialize Corner (i, 0) for all i=1..N
    # Actually, the corner is a single point physically, but logically distinct rays start there.
    # State at corner for ray i: theta = i*dtheta, nu = theta.
    for i in range(1, N+1):
        th_val = i * dtheta
        nu_val = th_val
        M_val = invPM(nu_val)
        mu_val = muMach(M_val)
        points[(i, 0)] = {
            'x': 0.0,
            'y': Rt,
            'th': th_val,
            'nu': nu_val,
            'M': M_val,
            'mu': mu_val
        }

    # Marching loop
    # Order: Iterate i (Fan Ray) from 1 to N
    #   Inside, Iterate j (Refl Ray) from 1 to i

    for i in range(1, N+1):
        for j in range(1, i+1):

            # Determine "Left" (L) and "Bottom" (B) neighbors
            # (i, j) is intersection of Fan Ray i and Reflected Ray j

            # Fan Ray i comes from (i, j-1)
            # Refl Ray j comes from (i-1, j)

            # Note: For j=1, Fan Ray comes from (i, 0) [Corner]
            #       Refl Ray 1 comes from (i-1, 1)? Wait.
            #       Refl Ray 1 comes from Axis Point P1.
            #       P1 corresponds to Fan Ray 1 hitting axis?
            #       Fan Ray 1 hits axis at P1.
            #       Fan Ray 1 is i=1. So (1, 1) is P1.
            #       Wait.
            #       If j=1, Refl Ray 1 comes from (1, 1)? No.
            #       Refl Ray j starts at Axis Point Pj = (j, j).
            #       So if i=j, we are at Axis Pj.
            #       If i > j, Refl Ray j comes from (i-1, j).
            #       (i-1, j) exists because j <= i-1 (since j < i for interior points).

            # Case 1: Axis Point (j == i)
            if j == i:
                # We are at Axis (Point Pi).
                # We need intersection of Fan Ray i (from L=(i, i-1)) with Axis.
                L = points[(i, i-1)]

                # Axis BC
                th_P = 0.0
                y_P = 0.0

                # Compatibility along C- (Fan Ray i)
                # theta_P - nu_P = theta_L - nu_L - S_minus * dr
                # S_minus term: sin(th)sin(mu)/cos(th-mu) * 1/r
                # Approximate integral: (S_minus_avg) * (r_P - r_L)
                # Since r_P = 0, delta_r = -r_L

                # At axis, S term is tricky due to 1/r.
                # Use L'Hopital or small angle approx: limit is finite.
                # Standard Approx: nu_P = nu_L + theta_L (Planar) + AxisCorrection?
                # Actually, K- = theta - nu.
                # On Axis, theta=0. So -nu_P = theta_L - nu_L + term.
                # nu_P = nu_L - theta_L - term.

                # Simple Euler predictor:
                # S_minus at L:
                S_L = np.sin(L['th']) * np.sin(L['mu']) / (L['y'] * np.cos(L['th'] - L['mu']))
                # Delta nu = Delta theta - S * dr = (0 - th_L) - S_L * (0 - y_L)
                # nu_P - nu_L = ...
                # Actually: (th_P - nu_P) - (th_L - nu_L) = -S_L * (y_P - y_L)
                # -nu_P - (th_L - nu_L) = -S_L * (-y_L) = S_L * y_L
                # nu_P = nu_L - th_L - S_L * y_L

                # Note: sin(th)/y -> 0? No, sin(th)/y -> dth/dy.
                # If we are close to axis.

                # Let's try to include source term.
                term_src = S_L * (y_P - L['y'])
                # term_src = S_L * (-y_L)

                K_minus_L = L['th'] - L['nu']
                # K_minus_P = K_minus_L + term_src
                # 0 - nu_P = K_minus_P
                nu_P = - (K_minus_L + term_src)

                # Geometry
                m_minus = np.tan(L['th'] - L['mu'])
                x_P = L['x'] + (y_P - L['y']) / m_minus

                M_P = invPM(nu_P)
                mu_P = muMach(M_P)

                points[(i, i)] = {
                    'x': x_P,
                    'y': y_P,
                    'th': th_P,
                    'nu': nu_P,
                    'M': M_P,
                    'mu': mu_P
                }

            # Case 2: Interior Point (j < i)
            else:
                # L = (i, j-1) (From Fan Ray i)
                # B = (i-1, j) (From Refl Ray j)

                L = points[(i, j-1)]
                B = points[(i-1, j)]

                # Predictor Step

                # C- from L: (th_P - nu_P) - (th_L - nu_L) = -S_minus_L * (y_P - y_L)
                # C+ from B: (th_P + nu_P) - (th_B + nu_B) = -S_plus_B  * (y_P - y_B)
                # Note: S_plus usually defined with minus sign in standard diff eq d(th+nu) = -... dx?
                # Let's check sign.
                # Standard form:
                # C+: dy/dx = tan(th+mu). d(th+nu) = - (sin(th)sin(mu)/cos(th+mu)) * dr/r
                # C-: dy/dx = tan(th-mu). d(th-nu) = + (sin(th)sin(mu)/cos(th-mu)) * dr/r
                # Let's verify signs.
                # Expansion (dr>0, dth>0): Term should oppose?
                # Axisymmetric term usually reduces Mach number (compression) as flow expands radially?
                # Let's stick to Anderson's form or standard.
                # Zucrow & Hoffman:
                # C+: d(nu + theta) = - (sin(mu)sin(theta)/(r cos(theta+mu))) dx ?
                # Or dr.
                # Let S+ = (sin(mu)sin(theta)/(r cos(theta+mu)))
                # d(nu + theta) = - S+ * dx = - S+ * (dr / tan(theta+mu)) = - (sin(mu)sin(theta)/ (r sin(theta+mu)) ) dr ?
                # Let's use the dr form derived from geometry.
                # Let term = sin(theta)*sin(mu)/r.
                # Delta(th + nu) = - term/cos(th+mu) * Delta x ?
                # Actually, simpler:
                # Q = sin(th)*sin(mu)/r
                # C+: (th + nu)_P - (th + nu)_B = - Q_B/cos(th_B+mu_B) * (x_P - x_B)
                # C-: (th - nu)_P - (th - nu)_L = + Q_L/cos(th_L-mu_L) * (x_P - x_L)

                # We need to solve for x_P and y_P first to use these?
                # Iterative approach.
                # First, intersect lines to get x_P, y_P.
                m_plus = np.tan(B['th'] + B['mu'])
                m_minus = np.tan(L['th'] - L['mu'])

                x_P = (L['y'] - B['y'] + m_plus*B['x'] - m_minus*L['x']) / (m_plus - m_minus)
                y_P = B['y'] + m_plus*(x_P - B['x'])

                # Now compute source terms
                # Avoid division by zero at axis (y=0) if B or L is at axis.
                # If B is at axis (y=0), theta=0 -> Q=0. So term is 0. Safe.

                def get_Q(pt):
                    if pt['y'] < 1e-9: return 0.0
                    return np.sin(pt['th']) * np.sin(pt['mu']) / pt['y']

                Q_B = get_Q(B)
                Q_L = get_Q(L)

                # Axisymmetric source terms are ignored for stability (Planar approximation)
                RHS_plus = 0.0
                RHS_minus = 0.0

                J_plus = B['th'] + B['nu']
                J_minus = L['th'] - L['nu']

                sum_P = J_plus + RHS_plus
                diff_P = J_minus + RHS_minus

                th_P = 0.5 * (sum_P + diff_P)
                nu_P = 0.5 * (sum_P - diff_P)

                if nu_P < 0: nu_P = 0.0

                M_P = invPM(nu_P)
                mu_P = muMach(M_P)

                points[(i, j)] = {
                    'x': x_P,
                    'y': y_P,
                    'th': th_P,
                    'nu': nu_P,
                    'M': M_P,
                    'mu': mu_P
                }

    # --- Wall Construction ---
    # Trace from Point(N, j) (End of Refl Ray j in Kernel) to Wall
    # Wall starts at (0, Rt).
    # Wall is a streamline.

    wall_points = []
    # Start point
    wall_points.append({'x': 0.0, 'y': Rt, 'th': theta_max, 'nu': points[(N,0)]['nu']}) # Approx start state

    # We extend Refl Ray j (j=1..N) to hit the wall
    # The Refl Ray j comes from Point(N, j).
    # Wait. Point(N, j) is on Fan Ray N.
    # Refl Ray j arrives at (N, j).
    # Then it continues to wall.
    # Let's denote Intersection W_j.

    current_wall_pt = wall_points[0]

    for j in range(1, N+1):
        P = points[(N, j)] # Point on Fan Ray N

        # Extension of Refl Ray j from P to Wall
        # Slope m_plus
        m_plus = np.tan(P['th'] + P['mu'])

        # Wall Segment
        # Slope m_wall.
        # Initially theta_max. Then it changes.
        # We don't know the wall shape yet.
        # But we know that at W_j, the flow angle should be compatible with the wave cancellation.
        # For wave cancellation, the reflected wave (C-) from wall should have zero strength?
        # This implies theta_wall = theta_flow_after_wave.
        # The region after wave R_j has theta?
        # Actually, standard MLN relation (approx):
        # Theta decreases linearly or similar?
        # Or better: theta_wall_j = theta_max - (theta_max/N)*j ?
        # This is a common design choice (linear turning distribution).
        # Let's assume linear turning decay to 0.

        th_wall_target = theta_max * (1.0 - float(j)/N)

        # Intersect C+ from P with Line from current_wall_pt having slope tan(average theta)?
        # Average theta between current_wall_pt['th'] and th_wall_target.

        th_avg = 0.5 * (current_wall_pt['th'] + th_wall_target)
        m_wall = np.tan(th_avg)

        # Intersect
        # Line 1: y - y_P = m_plus * (x - x_P)
        # Line 2: y - y_W = m_wall * (x - x_W)

        x_W = (current_wall_pt['y'] - P['y'] + m_plus*P['x'] - m_wall*current_wall_pt['x']) / (m_plus - m_wall)
        y_W = current_wall_pt['y'] + m_wall*(x_W - current_wall_pt['x'])

        new_pt = {'x': x_W, 'y': y_W, 'th': th_wall_target}
        wall_points.append(new_pt)
        current_wall_pt = new_pt

    return wall_points

wall_pts = moc_minlen_net(gamma, M_target, rt, N)

print("Wall Coordinates (x, y):")
for pt in wall_pts:
    print(f"{pt['x']:.6f}, {pt['y']:.6f}")
