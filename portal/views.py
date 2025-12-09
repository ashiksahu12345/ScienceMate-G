from django.shortcuts import render, redirect
import math

from .forms import ChatForm, QuestionForm
from .models import Question


# ===================== HOME: Subject + Topic buttons =====================

def home(request):
    """
    Main 'Science Mate' page.
    - Subject tabs (Physics / Maths / Chemistry)
    - Below: ONLY selected subject's topic buttons
    """

    subject = request.GET.get("subject", "physics").lower()
    if subject not in ["physics", "maths", "chemistry"]:
        subject = "physics"

    # ---------- Physics: full syllabus topic buttons ----------
    physics_topics = [
        {"label": "Units, Dimensions & Errors", "url": "physics_units"},
        {"label": "Motion (Kinematics)",           "url": "physics_motion"},
        {"label": "Newton’s Laws & Forces",        "url": "physics_forces"},
        {"label": "Work, Energy & Power",          "url": "physics_work_energy"},
        {"label": "Momentum & Collisions",         "url": "physics_momentum"},
        {"label": "Circular Motion",               "url": "physics_circular"},
        {"label": "Projectile Motion",             "url": "physics_projectile"},
        {"label": "Gravitation",                   "url": "physics_gravity"},
        {"label": "Fluid Mechanics",               "url": "physics_fluids"},
        {"label": "Heat & Thermodynamics",         "url": "physics_heat"},
        {"label": "Waves & Sound",                 "url": "physics_waves"},
        {"label": "SHM & Oscillations",            "url": "physics_shm"},
        {"label": "Ray Optics",                    "url": "physics_ray_optics"},
        {"label": "Wave Optics",                   "url": "physics_wave_optics"},
        {"label": "Electrostatics",                "url": "physics_electrostatics"},
        {"label": "Current Electricity (DC)",      "url": "physics_electricity"},
        {"label": "Magnetism & EMI",               "url": "physics_magnetism"},
        {"label": "AC Circuits",                   "url": "physics_ac"},
        {"label": "Electromagnetic Waves",         "url": "physics_em_waves"},
        {"label": "Modern Physics",                "url": "physics_modern"},
    ]

    # ---------- Maths topics ----------
        # ---------- Maths topics (full syllabus buttons) ----------
    maths_topics = [
        # Basic number system & algebra (Class 10 style)
        {"label": "Real Numbers",                         "url": "maths_real_numbers"},
        {"label": "Polynomials",                          "url": "maths_polynomials"},
        {"label": "Linear Equations in Two Variables",    "url": "maths_linear_equations"},
        {"label": "Quadratic Equations",                  "url": "maths_quadratic_equations"},
        {"label": "Arithmetic Progressions (AP)",         "url": "maths_ap"},

        {"label": "Coordinate Geometry (Basics)",         "url": "maths_coordinate_basics"},
        {"label": "Trigonometry – Basics",                "url": "maths_trig_basic"},
        {"label": "Trigonometry – Applications",          "url": "maths_trig_applications"},
        {"label": "Circles & Tangents",                   "url": "maths_circles"},
        {"label": "Areas Related to Circles",             "url": "maths_area_circles"},
        {"label": "Surface Area & Volume",                "url": "maths_surface_volume"},

        # Class 11–12 style topics
        {"label": "Relations & Functions",                "url": "maths_relations_functions"},
        {"label": "Matrices",                             "url": "maths_matrices"},
        {"label": "Determinants",                         "url": "maths_determinants"},
        {"label": "Limits & Continuity",                  "url": "maths_limits_continuity"},
        {"label": "Differentiation",                      "url": "maths_differentiation"},
        {"label": "Applications of Derivatives",          "url": "maths_app_derivatives"},
        {"label": "Integration",                          "url": "maths_integration"},
        {"label": "Applications of Integrals",            "url": "maths_app_integrals"},
        {"label": "Differential Equations",               "url": "maths_diff_eq"},
        {"label": "Vector Algebra",                       "url": "maths_vectors"},
        {"label": "3D Geometry",                          "url": "maths_3d_geometry"},
        {"label": "Linear Programming",                   "url": "maths_linear_programming"},

        # Existing calculator-style pages (already coded by you)
        {"label": "Arithmetic (Calculator)",              "url": "maths_arithmetic"},
        {"label": "Algebra (Linear & Quadratic Solver)",  "url": "maths_algebra"},
        {"label": "Geometry (Area Calculator)",           "url": "maths_geometry"},
        {"label": "Statistics (Mean/Variance/SD)",        "url": "maths_statistics"},
        {"label": "Probability (Basic)",                  "url": "maths_probability"},
        {"label": "Sphere (Surface Area & Volume)",       "url": "maths_sphere"},
        {"label": "Midpoint & Distance (Coordinate Geo)", "url": "maths_midpoint"},
    ]


    # ---------- Chemistry topics ----------
    chemistry_topics = [
        {"label": "Atomic Structure",   "url": "chem_atomic_structure"},
        {"label": "Chemical Bonding",   "url": "chem_chemical_bonding"},
        {"label": "Solutions",          "url": "chem_solutions"},
        {"label": "Thermodynamics",     "url": "chem_thermodynamics"},
    ]

    topics_map = {
        "physics": physics_topics,
        "maths": maths_topics,
        "chemistry": chemistry_topics,
    }

    context = {
        "active_subject": subject,
        "selected_subject": subject,
        "topics": topics_map[subject],
    }
    return render(request, "portal/home.html", context)


# ===================== Generic topic placeholder (for now) =====================

def topic_placeholder(request, subject, topic, description):
    """
    Generic page with title + short description.
    Later parts me yahi topics interactive calculators ban jayenge.
    """
    context = {
        "active_subject": subject,
        "subject_title": subject.capitalize(),
        "topic_title": topic,
        "description": description,
    }
    return render(request, "portal/topic_placeholder.html", context)


# ===================== PHYSICS — Full syllabus (placeholders) =====================

def physics_units(request):
    """Units, dimensions & error calculators page."""
    return render(request, "portal/physics_units.html")


# ===================== PHYSICS — INTERACTIVE TOPICS (PART 2) =====================

def physics_motion(request):
    """
    Motion (Kinematics):
    - v = u + a t
    - s = u t + ½ a t²
    - v² = u² + 2 a s
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "v_from_uat":
                u = float(request.POST["u"])
                a = float(request.POST["a"])
                t = float(request.POST["t"])
                v = u + a * t
                result = f"Final velocity v = {v:.4f} m/s"

            elif calc == "s_from_uat":
                u = float(request.POST["u"])
                a = float(request.POST["a"])
                t = float(request.POST["t"])
                s = u * t + 0.5 * a * t * t
                result = f"Displacement s = {s:.4f} m"

            elif calc == "v_from_uas":
                u = float(request.POST["u"])
                a = float(request.POST["a"])
                s = float(request.POST["s"])
                v2 = u * u + 2 * a * s
                if v2 < 0:
                    error = "v² is negative. Check values of a and s."
                else:
                    v = math.sqrt(v2)
                    result = f"Final velocity v = {v:.4f} m/s"
        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_motion.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_forces(request):
    """
    Newton’s Laws & Forces:
    - F = m a
    - Weight W = m g
    - Friction F = μ N
    """
    result = None
    error = None
    g = 9.8

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "f_ma":
                m = float(request.POST["m"])
                a = float(request.POST["a"])
                F = m * a
                result = f"Force F = {F:.4f} N"

            elif calc == "weight":
                m = float(request.POST["m"])
                W = m * g
                result = f"Weight W = {W:.4f} N"

            elif calc == "friction":
                mu = float(request.POST["mu"])
                N = float(request.POST["N"])
                F = mu * N
                result = f"Frictional force F = {F:.4f} N"
        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_forces.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_work_energy(request):
    """
    Work, Energy & Power:
    - Work: W = F d cosθ
    - Power: P = W / t
    - Kinetic Energy: KE = ½ m v²
    - Potential Energy: PE = m g h
    """
    result = None
    error = None
    g = 9.8

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # Work
            if calc == "work":
                F = float(request.POST["F"])
                d = float(request.POST["d"])
                theta = float(request.POST["theta"])
                W = F * d * math.cos(math.radians(theta))
                result = f"Work done W = {W:.4f} J"

            # Power
            elif calc == "power":
                W = float(request.POST["W"])
                t = float(request.POST["t"])
                if t == 0:
                    error = "Time t cannot be zero."
                else:
                    P = W / t
                    result = f"Power P = {P:.4f} W"

            # KE
            elif calc == "ke":
                m = float(request.POST["m"])
                v = float(request.POST["v"])
                KE = 0.5 * m * v * v
                result = f"Kinetic Energy KE = {KE:.4f} J"

            # PE
            elif calc == "pe":
                m = float(request.POST["m"])
                h = float(request.POST["h"])
                PE = m * g * h
                result = f"Potential Energy PE = {PE:.4f} J"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_work_energy.html",
        {"active_subject": "physics", "result": result, "error": error},
    )



# ===================== PHYSICS — PART 3 INTERACTIVE TOPICS =====================

def physics_momentum(request):
    """
    Momentum & Collisions:
    - p = m v
    - Impulse J = F Δt
    - Impulse J = m Δv
    (Abhi simple calculators; collisions ka main logic theory form me rahega.)
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # p = m v
            if calc == "p_mv":
                m = float(request.POST["m"])
                v = float(request.POST["v"])
                p = m * v
                result = f"Momentum p = {p:.4f} kg·m/s"

            # J = F Δt
            elif calc == "impulse_ft":
                F = float(request.POST["F"])
                dt = float(request.POST["dt"])
                J = F * dt
                result = f"Impulse J = {J:.44f} N·s"

            # J = m Δv
            elif calc == "impulse_mdv":
                m = float(request.POST["m"])
                dv = float(request.POST["dv"])
                J = m * dv
                result = f"Impulse (Δp) = {J:.4f} kg·m/s"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_momentum.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_gravity(request):
    """
    Gravitation:
    - Universal law: F = G m₁ m₂ / r²
    - Surface g:    g = G M / R²
    """
    result = None
    error = None
    G = 6.67e-11  # N·m²/kg²

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # F = G m1 m2 / r^2
            if calc == "grav_force":
                m1 = float(request.POST["m1"])
                m2 = float(request.POST["m2"])
                r = float(request.POST["r"])
                if r == 0:
                    error = "Distance r cannot be zero."
                else:
                    F = G * m1 * m2 / (r * r)
                    result = f"Gravitational force F = {F:.4e} N"

            # g = G M / R^2
            elif calc == "g_surface":
                M = float(request.POST["M"])
                R = float(request.POST["R"])
                if R == 0:
                    error = "Radius R cannot be zero."
                else:
                    g = G * M / (R * R)
                    result = f"Acceleration due to gravity g = {g:.4f} m/s²"

        except (KeyError, ValueError, ZeroDivisionError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_gravity.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_circular(request):
    """
    Circular Motion:
    - v = ω r
    - a_c = v² / r
    - a_c = ω² r
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # v = ω r
            if calc == "v_from_omega_r":
                omega = float(request.POST["omega"])
                r = float(request.POST["r"])
                v = omega * r
                result = f"Linear speed v = {v:.4f} m/s"

            # a_c = v² / r
            elif calc == "ac_from_vr":
                v = float(request.POST["v"])
                r = float(request.POST["r"])
                if r == 0:
                    error = "Radius r cannot be zero."
                else:
                    ac = v * v / r
                    result = f"Centripetal acceleration a_c = {ac:.4f} m/s²"

            # a_c = ω² r
            elif calc == "ac_from_omega_r":
                omega = float(request.POST["omega"])
                r = float(request.POST["r"])
                ac = omega * omega * r
                result = f"Centripetal acceleration a_c = {ac:.4f} m/s²"

        except (KeyError, ValueError, ZeroDivisionError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_circular.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_projectile(request):
    """
    Projectile Motion on level ground:
    Inputs: u (initial speed), θ (launch angle in degrees)
    Outputs:
      Time of flight: T = 2 u sinθ / g
      Range:         R = u² sin(2θ) / g
      Max height:    H = u² sin²θ / (2 g)
    """
    result = None
    error = None
    g = 9.8

    if request.method == "POST":
        try:
            u = float(request.POST["u"])
            angle = float(request.POST["angle"])
            rad = math.radians(angle)

            T = 2 * u * math.sin(rad) / g
            R = (u * u * math.sin(2 * rad)) / g
            H = (u * u * (math.sin(rad) ** 2)) / (2 * g)

            result = (
                f"Time of flight T = {T:.4f} s, "
                f"Range R = {R:.4f} m, "
                f"Maximum height H = {H:.4f} m"
            )
        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_projectile.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


# ===================== PHYSICS — PART 4 INTERACTIVE TOPICS =====================

def physics_fluids(request):
    """
    Fluid Mechanics:
    - Pressure in a liquid: P = ρ g h
    - Buoyant force: F_b = ρ g V
    - Continuity equation: A₁ v₁ = A₂ v₂  → v₂ = (A₁ v₁) / A₂
    """
    result = None
    error = None
    g = 9.8

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "pressure":
                rho = float(request.POST["rho"])
                h = float(request.POST["h"])
                P = rho * g * h
                result = f"Pressure P = {P:.4f} Pa"

            elif calc == "buoyant":
                rho = float(request.POST["rho"])
                V = float(request.POST["V"])
                Fb = rho * g * V
                result = f"Buoyant force Fᵦ = {Fb:.4f} N"

            elif calc == "continuity":
                A1 = float(request.POST["A1"])
                v1 = float(request.POST["v1"])
                A2 = float(request.POST["A2"])
                if A2 == 0:
                    error = "A₂ cannot be zero."
                else:
                    v2 = (A1 * v1) / A2
                    result = f"Final speed v₂ = {v2:.4f} m/s"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_fluids.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_heat(request):
    """
    Heat & Thermodynamics (Physics):
    - Heat: Q = m c ΔT
    - Linear expansion: ΔL = α L ΔT
    - First law: ΔU = Q - W
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # Q = m c ΔT
            if calc == "heat":
                m = float(request.POST["m"])
                c = float(request.POST["c"])
                dT = float(request.POST["dT"])
                Q = m * c * dT
                result = f"Heat absorbed/released Q = {Q:.4f} J"

            # ΔL = α L ΔT
            elif calc == "expansion":
                alpha = float(request.POST["alpha"])
                L = float(request.POST["L"])
                dT = float(request.POST["dT"])
                dL = alpha * L * dT
                result = f"Change in length ΔL = {dL:.6f} m"

            # ΔU = Q - W
            elif calc == "delta_u":
                Q = float(request.POST["Q"])
                W = float(request.POST["W"])
                dU = Q - W
                result = f"Change in internal energy ΔU = {dU:.4f} J"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_heat.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_waves(request):
    """
    Waves & Sound:
    - Wave speed: v = f λ
    - Time period: T = 1 / f
    - Speed of sound in air (approx): v ≈ 331 + 0.6 T(°C)
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # v = f λ
            if calc == "v_fl":
                f = float(request.POST["f"])
                lam = float(request.POST["lam"])
                v = f * lam
                result = f"Wave speed v = {v:.4f} m/s"

            # T = 1 / f
            elif calc == "t_from_f":
                f = float(request.POST["f"])
                if f == 0:
                    error = "Frequency f cannot be zero."
                else:
                    T = 1.0 / f
                    result = f"Time period T = {T:.6f} s"

            # v ≈ 331 + 0.6 T
            elif calc == "sound_speed":
                temp = float(request.POST["temp"])
                v = 331.0 + 0.6 * temp
                result = f"Approx. speed of sound v ≈ {v:.2f} m/s"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_waves.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_shm(request):
    """
    SHM & Oscillations:
    - Spring–mass system: T = 2π √(m/k)
    - Simple pendulum:    T = 2π √(L/g)
    - Max speed in SHM:   v_max = ω A = (2π/T) A
    """
    result = None
    error = None
    g = 9.8

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # Spring–mass T
            if calc == "spring_T":
                m = float(request.POST["m"])
                k = float(request.POST["k"])
                if k <= 0:
                    error = "Spring constant k must be > 0."
                else:
                    T = 2 * math.pi * math.sqrt(m / k)
                    result = f"Time period (spring–mass) T = {T:.4f} s"

            # Pendulum T
            elif calc == "pendulum_T":
                L = float(request.POST["L"])
                if L <= 0:
                    error = "Length L must be > 0."
                else:
                    T = 2 * math.pi * math.sqrt(L / g)
                    result = f"Time period (pendulum) T = {T:.4f} s"

            # v_max = (2π/T) A
            elif calc == "vmax":
                A = float(request.POST["A"])
                T = float(request.POST["T"])
                if T <= 0:
                    error = "Time period T must be > 0."
                else:
                    omega = 2 * math.pi / T
                    vmax = omega * A
                    result = f"Maximum speed v_max = {vmax:.4f} (units/s)"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_shm.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_ray_optics(request):
    """Ray optics calculators page."""
    return render(request, "portal/physics_ray_optics.html")



def physics_wave_optics(request):
    """Wave optics calculators page."""
    return render(request, "portal/physics_wave_optics.html")


# ===================== PHYSICS — PART 5 INTERACTIVE TOPICS =====================

def physics_electrostatics(request):
    """
    Electrostatics:
    - Coulomb force:      F = k q₁ q₂ / r²
    - Electric field:     E = F / q
    - Potential energy:   U = k q₁ q₂ / r
    Units here:
      q in Coulomb, r in metre, k ≈ 9×10⁹ N·m²/C²
    """
    result = None
    error = None
    k = 9e9

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "coulomb":
                q1 = float(request.POST["q1"])
                q2 = float(request.POST["q2"])
                r = float(request.POST["r"])
                if r == 0:
                    error = "Distance r cannot be zero."
                else:
                    F = k * q1 * q2 / (r * r)
                    result = f"Coulomb force F = {F:.4e} N"

            elif calc == "field":
                F = float(request.POST["F"])
                q = float(request.POST["q"])
                if q == 0:
                    error = "Charge q cannot be zero."
                else:
                    E = F / q
                    result = f"Electric field E = {E:.4e} N/C"

            elif calc == "energy":
                q1 = float(request.POST["q1"])
                q2 = float(request.POST["q2"])
                r = float(request.POST["r"])
                if r == 0:
                    error = "Distance r cannot be zero."
                else:
                    U = k * q1 * q2 / r
                    result = f"Potential energy U = {U:.4e} J"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_electrostatics.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_electricity(request):
    """
    Current Electricity (DC):
    - Ohm’s law:        V = I R
    - Power:            P = V I
    - R_series:         R_eq = R₁ + R₂ + R₃
    - R_parallel:       1/R_eq = 1/R₁ + 1/R₂ + 1/R₃
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # Ohm's law
            if calc == "ohm":
                I = float(request.POST["I"])
                R = float(request.POST["R"])
                V = I * R
                result = f"Voltage V = {V:.4f} V"

            # Power
            elif calc == "power":
                V = float(request.POST["V"])
                I = float(request.POST["I"])
                P = V * I
                result = f"Power P = {P:.4f} W"

            # R series
            elif calc == "r_series":
                R1 = float(request.POST["R1"])
                R2 = float(request.POST["R2"])
                R3 = float(request.POST["R3"])
                Req = R1 + R2 + R3
                result = f"Equivalent resistance (series) R_eq = {Req:.4f} Ω"

            # R parallel
            elif calc == "r_parallel":
                R1 = float(request.POST["R1"])
                R2 = float(request.POST["R2"])
                R3 = float(request.POST["R3"])
                if R1 == 0 or R2 == 0 or R3 == 0:
                    error = "For this simple formula, none of R₁, R₂, R₃ should be zero."
                else:
                    Req_inv = 1.0 / R1 + 1.0 / R2 + 1.0 / R3
                    Req = 1.0 / Req_inv
                    result = f"Equivalent resistance (parallel) R_eq = {Req:.4f} Ω"

        except (KeyError, ValueError, ZeroDivisionError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_electricity.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_magnetism(request):
    """
    Magnetism & Electromagnetic Induction:
    - Force on charge:       F = q v B sinθ
    - Force on conductor:    F = B I L sinθ
    - Induced emf:           |ε| = N ΔΦ / Δt
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # q v B sinθ
            if calc == "force_charge":
                q = float(request.POST["q"])
                v = float(request.POST["v"])
                B = float(request.POST["B"])
                theta = float(request.POST["theta"])
                F = q * v * B * math.sin(math.radians(theta))
                result = f"Magnetic force on charge F = {F:.4e} N"

            # B I L sinθ
            elif calc == "force_wire":
                I = float(request.POST["I"])
                L = float(request.POST["L"])
                B = float(request.POST["B"])
                theta = float(request.POST["theta"])
                F = B * I * L * math.sin(math.radians(theta))
                result = f"Force on wire F = {F:.4e} N"

            # ε = N ΔΦ / Δt
            elif calc == "emf":
                N = float(request.POST["N"])
                dphi = float(request.POST["dphi"])
                dt = float(request.POST["dt"])
                if dt == 0:
                    error = "Time interval Δt cannot be zero."
                else:
                    emf = N * dphi / dt
                    result = f"Induced emf (magnitude) ε = {emf:.4e} V"

        except (KeyError, ValueError, ZeroDivisionError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_magnetism.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_ac(request):
    """
    AC Circuits:
    - V_rms = V_max / √2 and I_rms = I_max / √2
    - Power in AC: P = V_rms I_rms cosφ
    - Reactances: X_L = 2π f L,  X_C = 1 / (2π f C)
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # Vrms & Irms
            if calc == "rms":
                Vmax = float(request.POST["Vmax"])
                Imax = float(request.POST["Imax"])
                Vrms = Vmax / math.sqrt(2)
                Irms = Imax / math.sqrt(2)
                result = f"V_rms = {Vrms:.4f} V,  I_rms = {Irms:.4f} A"

            # AC power
            elif calc == "power_ac":
                Vr = float(request.POST["Vr"])
                Ir = float(request.POST["Ir"])
                pf = float(request.POST["pf"])
                P = Vr * Ir * pf
                result = f"Average power P = {P:.4f} W"

            # Reactances
            elif calc == "reactance":
                f = float(request.POST["f"])
                L = float(request.POST["L"])
                C = float(request.POST["C"])
                if f <= 0:
                    error = "Frequency f must be > 0."
                else:
                    XL = 2 * math.pi * f * L
                    XC = 0.0
                    if C > 0:
                        XC = 1.0 / (2 * math.pi * f * C)
                    result = (
                        f"Inductive reactance X_L = {XL:.4f} Ω, "
                        f"Capacitive reactance X_C = {XC:.4f} Ω"
                    )

        except (KeyError, ValueError, ZeroDivisionError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_ac.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_em_waves(request):
    """
    Electromagnetic Waves:
    - c = λ f  → f = c / λ or λ = c / f
    - Energy of photon: E = h f
    (c ≈ 3×10⁸ m/s, h ≈ 6.626×10⁻³⁴ J·s)
    """
    result = None
    error = None
    c = 3e8
    h = 6.626e-34

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "freq_from_lambda":
                lam = float(request.POST["lam"])
                if lam <= 0:
                    error = "Wavelength λ must be > 0."
                else:
                    f = c / lam
                    result = f"Frequency f = {f:.4e} Hz"

            elif calc == "lambda_from_freq":
                f = float(request.POST["f"])
                if f <= 0:
                    error = "Frequency f must be > 0."
                else:
                    lam = c / f
                    result = f"Wavelength λ = {lam:.4e} m"

            elif calc == "photon_energy":
                f = float(request.POST["f"])
                if f <= 0:
                    error = "Frequency f must be > 0."
                else:
                    E = h * f
                    result = f"Photon energy E = {E:.4e} J"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_em_waves.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_modern(request):
    """
    Modern Physics:
    - Mass–energy:         E = m c²
    - de Broglie λ:        λ = h / (m v)
    - Nuclear energy (simple): E = Δm c²
    """
    result = None
    error = None
    c = 3e8
    h = 6.626e-34

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # E = m c^2
            if calc == "mass_energy":
                m = float(request.POST["m"])
                E = m * c * c
                result = f"Energy E = {E:.4e} J"

            # λ = h / (m v)
            elif calc == "de_broglie":
                m = float(request.POST["m"])
                v = float(request.POST["v"])
                if m <= 0 or v <= 0:
                    error = "Mass and velocity must both be > 0."
                else:
                    lam = h / (m * v)
                    result = f"De Broglie wavelength λ = {lam:.4e} m"

            # Nuclear Δm c^2
            elif calc == "nuclear":
                dm = float(request.POST["dm"])
                E = dm * c * c
                result = f"Released energy E = {E:.4e} J"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_modern.html",
        {"active_subject": "physics", "result": result, "error": error},
    )



# ===================== MATHS — Interactive (as built earlier) =====================

def maths_arithmetic(request):
    """
    Arithmetic:
    - Percentage = (part / whole) × 100
    - Average of n numbers (comma separated)
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "percentage":
                part = float(request.POST["part"])
                whole = float(request.POST["whole"])
                if whole == 0:
                    error = "Whole cannot be zero."
                else:
                    pct = (part / whole) * 100
                    result = f"Percentage = {pct:.2f}%"

            elif calc == "average":
                raw = request.POST["numbers"]
                nums = [float(x) for x in raw.replace(" ", "").split(",") if x != ""]
                if not nums:
                    error = "Please enter at least one number."
                else:
                    avg = sum(nums) / len(nums)
                    result = f"Average = {avg:.4f}"
        except (ValueError, KeyError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_arithmetic.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_algebra(request):
    """
    Algebra:
    - Linear: ax + b = 0
    - Quadratic: ax² + bx + c = 0
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "linear":
                a = float(request.POST["a"])
                b = float(request.POST["b"])
                if a == 0:
                    error = "For linear equation, a cannot be 0."
                else:
                    x = -b / a
                    result = f"Solution: x = {x:.4f}"

            elif calc == "quadratic":
                a = float(request.POST["a"])
                b = float(request.POST["b"])
                c = float(request.POST["c"])
                if a == 0:
                    error = "This is not a quadratic equation (a = 0)."
                else:
                    D = b * b - 4 * a * c
                    if D > 0:
                        x1 = (-b + math.sqrt(D)) / (2 * a)
                        x2 = (-b - math.sqrt(D)) / (2 * a)
                        result = f"Two real roots: x₁ = {x1:.4f}, x₂ = {x2:.4f}"
                    elif D == 0:
                        x = -b / (2 * a)
                        result = f"One real root: x = {x:.4f}"
                    else:
                        error = "No real roots (discriminant < 0)."
        except (ValueError, KeyError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_algebra.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_geometry(request):
    """
    Geometry:
    - Area of rectangle: A = l × b
    - Area of triangle: A = ½ b h
    - Area of circle: A = π r²
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "rect":
                l = float(request.POST["l"])
                b = float(request.POST["b"])
                A = l * b
                result = f"Area of rectangle = {A:.4f} square units"

            elif calc == "tri":
                base = float(request.POST["b"])
                h = float(request.POST["h"])
                A = 0.5 * base * h
                result = f"Area of triangle = {A:.4f} square units"

            elif calc == "circle":
                r = float(request.POST["r"])
                A = math.pi * r * r
                result = f"Area of circle = {A:.4f} square units"
        except (ValueError, KeyError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_geometry.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_statistics(request):
    """
    Statistics:
    - Mean, variance, standard deviation of comma-separated numbers
    """
    result = None
    error = None

    if request.method == "POST":
        try:
            raw = request.POST["numbers"]
            nums = [float(x) for x in raw.replace(" ", "").split(",") if x != ""]
            if not nums:
                error = "Please enter at least one number."
            else:
                n = len(nums)
                mean = sum(nums) / n
                variance = sum((x - mean) ** 2 for x in nums) / n
                stddev = math.sqrt(variance)
                result = (
                    f"Mean = {mean:.4f}, Variance = {variance:.4f}, "
                    f"Standard Deviation = {stddev:.4f}"
                )
        except (ValueError, KeyError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_statistics.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_probability(request):
    """
    Probability:
    - P(A) = favourable / total
    - P(A') = 1 - P(A)
    """
    result = None
    error = None

    if request.method == "POST":
        try:
            fav = float(request.POST["fav"])
            total = float(request.POST["total"])
            if total <= 0:
                error = "Total outcomes must be > 0."
            elif fav < 0 or fav > total:
                error = "Favourable outcomes must be between 0 and total."
            else:
                P = fav / total
                Pc = 1 - P
                result = f"P(A) = {P:.4f}, P(A') = {Pc:.4f}"
        except (ValueError, KeyError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_probability.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_sphere(request):
    """
    Sphere:
    - Surface area: SA = 4πr²
    - Volume: V = 4/3 π r³
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            r = float(request.POST["r"])
            if r < 0:
                error = "Radius cannot be negative."
            else:
                if calc == "sa":
                    SA = 4 * math.pi * r * r
                    result = f"Surface area SA = {SA:.4f} square units"
                elif calc == "vol":
                    V = (4.0 / 3.0) * math.pi * r ** 3
                    result = f"Volume V = {V:.4f} cubic units"
        except (ValueError, KeyError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_sphere.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_midpoint(request):
    """
    Coordinate Geometry:
    - Midpoint: M = ((x1 + x2)/2, (y1 + y2)/2)
    - Distance: d = √((x2 - x1)² + (y2 - y1)²)
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            x1 = float(request.POST["x1"])
            y1 = float(request.POST["y1"])
            x2 = float(request.POST["x2"])
            y2 = float(request.POST["y2"])

            if calc == "mid":
                mx = (x1 + x2) / 2.0
                my = (y1 + y2) / 2.0
                result = f"Midpoint M = ({mx:.4f}, {my:.4f})"
            elif calc == "dist":
                d = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
                result = f"Distance d = {d:.4f}"
        except (ValueError, KeyError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_midpoint.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


# ===================== CHEMISTRY — Interactive =====================

def chem_atomic_structure(request):
    """
    Atomic Structure:
    - Bohr energy: Eₙ = -13.6 Z² / n²  (eV)
    - De Broglie wavelength: λ = h / (m v)
    """
    result = None
    error = None
    h = 6.626e-34  # Planck's constant (J·s)

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "bohr_energy":
                Z = float(request.POST["Z"])
                n = float(request.POST["n"])
                if n <= 0:
                    error = "n must be > 0."
                else:
                    E = -13.6 * (Z ** 2) / (n ** 2)
                    result = f"Energy of electron Eₙ = {E:.4f} eV"

            elif calc == "de_broglie":
                m = float(request.POST["m"])
                v = float(request.POST["v"])
                if m <= 0 or v <= 0:
                    error = "Mass and velocity must be > 0."
                else:
                    lamb = h / (m * v)
                    result = f"De Broglie wavelength λ = {lamb:.4e} m"
        except (ValueError, KeyError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/chem_atomic_structure.html",
        {"active_subject": "chemistry", "result": result, "error": error},
    )


def chem_chemical_bonding(request):
    """
    Chemical Bonding:
    - Approx % ionic character (Pauling relation):
      %IC ≈ 16(Δχ) + 3.5(Δχ)²
      where Δχ = |χA - χB|
    """
    result = None
    error = None

    if request.method == "POST":
        try:
            chiA = float(request.POST["chiA"])
            chiB = float(request.POST["chiB"])
            dx = abs(chiA - chiB)
            pct = 16 * dx + 3.5 * dx * dx
            result = f"Approx. % ionic character ≈ {pct:.2f}%"
        except (ValueError, KeyError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/chem_chemical_bonding.html",
        {"active_subject": "chemistry", "result": result, "error": error},
    )


def chem_solutions(request):
    """
    Solutions:
    - Molarity M = n / V
    - Mass percent = (mass solute / total mass) × 100
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "molarity":
                n = float(request.POST["n"])
                V = float(request.POST["V"])
                if V == 0:
                    error = "Volume V cannot be zero."
                else:
                    M = n / V
                    result = f"Molarity M = {M:.4f} mol/L"

            elif calc == "mass_percent":
                ms = float(request.POST["ms"])
                msolvent = float(request.POST["msolvent"])
                total = ms + msolvent
                if total == 0:
                    error = "Total mass cannot be zero."
                else:
                    pct = (ms / total) * 100
                    result = f"Mass percent = {pct:.2f}%"
        except (ValueError, KeyError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/chem_solutions.html",
        {"active_subject": "chemistry", "result": result, "error": error},
    )


def chem_thermodynamics(request):
    """
    Thermodynamics:
    - First law: ΔU = Q - W
    - Gibbs free energy: ΔG = ΔH - T ΔS
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "delta_u":
                Q = float(request.POST["Q"])
                W = float(request.POST["W"])
                dU = Q - W
                result = f"Change in internal energy ΔU = {dU:.4f} J"

            elif calc == "delta_g":
                dH = float(request.POST["dH"])
                T = float(request.POST["T"])
                dS = float(request.POST["dS"])
                dG = dH - T * dS
                result = f"Change in Gibbs free energy ΔG = {dG:.4f} J"
        except (ValueError, KeyError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/chem_thermodynamics.html",
        {"active_subject": "chemistry", "result": result, "error": error},
    )


# ===================== CHATBOT (same as earlier) =====================

def _generate_bot_reply(user_text: str) -> str:
    text = user_text.lower()

    if "momentum" in text:
        return "Momentum p = m × v and impulse J = Δp = F × Δt."
    if "sphere" in text:
        return "Surface area of a sphere: SA = 4πr². Volume = 4/3 πr³."
    if "midpoint" in text:
        return "Midpoint of two points: M = ((x₁ + x₂)/2 , (y₁ + y₂)/2)."
    if "probability" in text:
        return "Basic probability P = favourable outcomes / total outcomes."
    if "hello" in text or "hi" in text:
        return "Hi! I am Science Mate Bot. Ask me doubts from Physics, Maths or Chemistry."

    return (
        "Nice question! Right now I give small hints using formulas. "
        "For full solution, try topic pages or upload your question."
    )


def chatbot(request):
    chat_history = request.session.get("chat_history", [])

    if request.method == "POST":
        form = ChatForm(request.POST)
        if form.is_valid():
            user_msg = form.cleaned_data["message"]
            bot_msg = _generate_bot_reply(user_msg)

            chat_history.append({"sender": "user", "text": user_msg})
            chat_history.append({"sender": "bot", "text": bot_msg})
            chat_history = chat_history[-30:]
            request.session["chat_history"] = chat_history
            return redirect("chatbot")
    else:
        form = ChatForm()
        if not chat_history:
            chat_history = [
                {"sender": "bot", "text": "Welcome to Science Mate Chatbot! 👋"},
                {
                    "sender": "bot",
                    "text": "Ask me about formulas like momentum, sphere, midpoint, probability etc.",
                },
            ]
            request.session["chat_history"] = chat_history

    context = {
        "form": form,
        "chat_history": chat_history,
        "active_subject": None,
    }
    return render(request, "portal/chatbot.html", context)


# ===================== UPLOAD QUESTION (same as earlier) =====================

def upload_question(request):
    if request.method == "POST":
        form = QuestionForm(request.POST)
        if form.is_valid():
            form.save()
            return redirect("upload_question")
    else:
        form = QuestionForm()

    recent_questions = Question.objects.order_by("-created_at")[:10]

    context = {
        "form": form,
        "recent_questions": recent_questions,
        "active_subject": None,
    }
    return render(request, "portal/upload_question.html", context)


def topic_placeholder(request, subject, topic, description):
    """
    Generic page with title + short description.
    Later parts me yahi topics interactive calculators ban jayenge.
    """
    context = {
        "active_subject": subject,
        "subject_title": subject.capitalize(),
        "topic_title": topic,
        "description": description,
    }
    return render(request, "portal/topic_placeholder.html", context)


# ===================== MATHS — Full syllabus placeholder pages =====================

def maths_real_numbers(request):
    return topic_placeholder(
        request,
        "maths",
        "Real Numbers",
        "Properties of integers, rational and irrational numbers, Euclid’s division lemma, HCF and LCM, etc."
    )


def maths_polynomials(request):
    return topic_placeholder(
        request,
        "maths",
        "Polynomials",
        "Degree of polynomial, zeros of polynomials, relationship between zeros and coefficients, graphs, factorisation."
    )


def maths_linear_equations(request):
    return topic_placeholder(
        request,
        "maths",
        "Linear Equations in Two Variables",
        "System of linear equations, graphical and algebraic solutions, elimination and substitution methods."
    )


def maths_quadratic_equations(request):
    return topic_placeholder(
        request,
        "maths",
        "Quadratic Equations",
        "Standard form ax² + bx + c = 0, discriminant, nature of roots, factorisation and quadratic formula."
    )


def maths_ap(request):
    return topic_placeholder(
        request,
        "maths",
        "Arithmetic Progressions (AP)",
        "nth term and sum of n terms of an AP, applications in word problems."
    )


def maths_coordinate_basics(request):
    return topic_placeholder(
        request,
        "maths",
        "Coordinate Geometry (Basics)",
        "Cartesian plane, distance formula, section formula, area of triangle using coordinates."
    )


def maths_trig_basic(request):
    return topic_placeholder(
        request,
        "maths",
        "Trigonometry – Basics",
        "Trigonometric ratios, identities, complementary angles, and simple problems."
    )


def maths_trig_applications(request):
    return topic_placeholder(
        request,
        "maths",
        "Trigonometry – Applications",
        "Heights and distances, angle of elevation and depression, simple real-life problems."
    )


def maths_circles(request):
    return topic_placeholder(
        request,
        "maths",
        "Circles & Tangents",
        "Basic terms, tangent properties, number of tangents from a point, angle theorems."
    )


def maths_area_circles(request):
    return topic_placeholder(
        request,
        "maths",
        "Areas Related to Circles",
        "Area and perimeter of circle, sector and segment, combination of plane figures."
    )


def maths_surface_volume(request):
    return topic_placeholder(
        request,
        "maths",
        "Surface Area & Volume",
        "Surface area and volume of cuboid, cylinder, cone, sphere, frustum, and combined solids."
    )


def maths_relations_functions(request):
    return topic_placeholder(
        request,
        "maths",
        "Relations & Functions",
        "Ordered pairs, Cartesian product, relation, function, domain, co-domain and range."
    )


def maths_matrices(request):
    return topic_placeholder(
        request,
        "maths",
        "Matrices",
        "Types of matrices, operations, transpose, determinant for small order matrices."
    )


def maths_determinants(request):
    return topic_placeholder(
        request,
        "maths",
        "Determinants",
        "Determinants of order 2 and 3, properties, area of triangle, solving linear equations (Cramer's rule)."
    )


def maths_limits_continuity(request):
    return topic_placeholder(
        request,
        "maths",
        "Limits & Continuity",
        "Intuitive idea of limit, basic limit rules, continuity of simple functions."
    )


def maths_differentiation(request):
    return topic_placeholder(
        request,
        "maths",
        "Differentiation",
        "Derivative as rate of change, basic formulas, product, quotient and chain rule."
    )


def maths_app_derivatives(request):
    return topic_placeholder(
        request,
        "maths",
        "Applications of Derivatives",
        "Increasing/decreasing functions, maxima-minima of simple functions, tangents and normals."
    )


def maths_integration(request):
    return topic_placeholder(
        request,
        "maths",
        "Integration",
        "Indefinite integrals, basic formulas, substitution, and simple standard results."
    )


def maths_app_integrals(request):
    return topic_placeholder(
        request,
        "maths",
        "Applications of Integrals",
        "Area under curves and between curves using definite integrals (basic level)."
    )


def maths_diff_eq(request):
    return topic_placeholder(
        request,
        "maths",
        "Differential Equations",
        "Order and degree, formation of differential equation, basic methods of solution (variable separable, etc.)."
    )


def maths_vectors(request):
    return topic_placeholder(
        request,
        "maths",
        "Vector Algebra",
        "Basic concepts of vectors, addition, scalar multiplication, dot product and cross product (intro level)."
    )


def maths_3d_geometry(request):
    return topic_placeholder(
        request,
        "maths",
        "3D Geometry",
        "Direction cosines, equation of line and plane (very basic introduction)."
    )


def maths_linear_programming(request):
    return topic_placeholder(
        request,
        "maths",
        "Linear Programming",
        "Formulation of LPP, graphical method of solution for two variables, feasible region and optimal solution."
    )


