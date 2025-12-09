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
        {"label": "Units, Dimensions & Errors",      "url": "physics_units"},
        {"label": "Motion (Kinematics)",            "url": "physics_motion"},
        {"label": "Newton’s Laws & Forces",         "url": "physics_forces"},
        {"label": "Work, Energy & Power",           "url": "physics_work_energy"},
        {"label": "Momentum & Collisions",          "url": "physics_momentum"},
        {"label": "Circular Motion",                "url": "physics_circular"},
        {"label": "Projectile Motion",              "url": "physics_projectile"},

        {"label": "Gravitation",                    "url": "physics_gravity"},
        {"label": "Fluid Mechanics",                "url": "physics_fluids"},
        {"label": "Heat & Thermodynamics",          "url": "physics_heat"},
        {"label": "Waves & Sound",                  "url": "physics_waves"},
        {"label": "SHM & Oscillations",             "url": "physics_shm"},
        {"label": "Ray Optics",                     "url": "physics_ray_optics"},
        {"label": "Wave Optics",                    "url": "physics_wave_optics"},

        {"label": "Electrostatics",                 "url": "physics_electrostatics"},
        {"label": "Current Electricity (DC)",       "url": "physics_electricity"},
        {"label": "Magnetism & EMI",                "url": "physics_magnetism"},
        {"label": "AC Circuits",                    "url": "physics_ac"},
        {"label": "Electromagnetic Waves",          "url": "physics_em_waves"},
        {"label": "Modern Physics",                 "url": "physics_modern"},

        {"label": "Rotational Motion & Rigid Bodies", "url": "physics_rotational"},
        {"label": "Properties of Solids & Elasticity","url": "physics_solids"},
        {"label": "Kinetic Theory of Gases",          "url": "physics_kinetic_theory"},
        {"label": "Semiconductors & Electronics",     "url": "physics_semiconductors"},
        {"label": "Communication Systems",            "url": "physics_communication"},

        # ===== 10 NEW EXTRA TOPICS =====
        {"label": "Thermal Properties of Matter",      "url": "physics_thermal_properties"},
        {"label": "Electric Potential & Capacitance",  "url": "physics_capacitance"},
        {"label": "Magnetic Effects of Current",       "url": "physics_magnetic_effects"},
        {"label": "Electromagnetic Induction (Detailed)", "url": "physics_emi_detailed"},
        {"label": "AC Power & Resonance",              "url": "physics_ac_power"},
        {"label": "Laser Physics",                     "url": "physics_lasers"},
        {"label": "Optical Fiber & Waveguides",        "url": "physics_optical_fiber"},
        {"label": "Nuclear Physics & Radioactivity",   "url": "physics_nuclear"},
        {"label": "Error & Data Analysis",             "url": "physics_errors_data"},
        {"label": "Instruments & Measurements",        "url": "physics_instruments"},
    ]



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
    # ---------- Chemistry topics: full syllabus ----------
        # ---------- Chemistry topics: full syllabus ----------
    chemistry_topics = [
        # Existing interactive pages
        {"label": "Atomic Structure",                      "url": "chem_atomic_structure"},
        {"label": "Chemical Bonding & Molecular Structure","url": "chem_chemical_bonding"},
        {"label": "Solutions",                             "url": "chem_solutions"},
        {"label": "Thermodynamics",                        "url": "chem_thermodynamics"},

        # Physical & General Chemistry
        {"label": "Stoichiometry & Mole Concept",          "url": "chem_stoichiometry"},
        {"label": "States of Matter (Gas, Liquid & Solid)","url": "chem_states_matter"},
        {"label": "Periodic Table & Periodicity",          "url": "chem_periodic_table"},
        {"label": "Chemical & Ionic Equilibrium",          "url": "chem_equilibrium"},
        {"label": "Redox Reactions",                       "url": "chem_redox"},
        {"label": "Electrochemistry",                      "url": "chem_electrochemistry"},
        {"label": "Chemical Kinetics",                     "url": "chem_chemical_kinetics"},
        {"label": "Surface Chemistry",                     "url": "chem_surface_chemistry"},

        # Inorganic Chemistry
        {"label": "Hydrogen & its Compounds",              "url": "chem_hydrogen"},
        {"label": "s-Block Elements",                      "url": "chem_s_block"},
        {"label": "p-Block Elements",                      "url": "chem_p_block"},
        {"label": "d- & f-Block Elements",                 "url": "chem_d_f_block"},
        {"label": "Coordination Compounds",                "url": "chem_coordination"},
        {"label": "Metallurgy & Qualitative Analysis",     "url": "chem_metallurgy"},
        {"label": "Environmental Chemistry",               "url": "chem_environmental"},

        # Organic Chemistry (11th + 12th)
        {"label": "Basic Concepts of Organic & IUPAC",     "url": "chem_organic_basics"},
        {"label": "Hydrocarbons (Alkane, Alkene, Alkyne)", "url": "chem_hydrocarbons"},
        {"label": "Haloalkanes & Haloarenes",              "url": "chem_haloalkanes"},
        {"label": "Alcohols, Phenols & Ethers",            "url": "chem_alcohols"},
        {"label": "Aldehydes, Ketones & Carboxylic Acids", "url": "chem_aldehydes"},
        {"label": "Amines",                                "url": "chem_amines"},
        {"label": "Biomolecules",                          "url": "chem_biomolecules"},
        {"label": "Polymers",                              "url": "chem_polymers"},
        {"label": "Chemistry in Everyday Life",            "url": "chem_everyday_life"},
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

# ===================== PHYSICS — Extra syllabus topics (placeholders) =====================

def physics_thermal_properties(request):
    return topic_placeholder(
        request,
        "physics",
        "Thermal Properties of Matter",
        "Expansion of solids, liquids and gases, heat capacity, specific heat, calorimetry, "
        "Newton’s law of cooling and basic ideas of thermal conductivity."
    )


def physics_capacitance(request):
    return topic_placeholder(
        request,
        "physics",
        "Electric Potential & Capacitance",
        "Electric potential and potential difference, equipotential surfaces, capacitance of parallel plate capacitor, "
        "energy stored in capacitor and combination of capacitors."
    )


def physics_magnetic_effects(request):
    return topic_placeholder(
        request,
        "physics",
        "Magnetic Effects of Current",
        "Biot–Savart law, Ampere’s law (basic applications), magnetic field due to current in straight wire, loop and solenoid, "
        "force on a moving charge and current carrying conductor."
    )


def physics_emi_detailed(request):
    return topic_placeholder(
        request,
        "physics",
        "Electromagnetic Induction (Detailed)",
        "Faraday’s law, Lenz’s law, induced emf and current, self and mutual inductance, "
        "induced emf in rotating coil and simple problems on EM induction."
    )


def physics_ac_power(request):
    return topic_placeholder(
        request,
        "physics",
        "AC Power & Resonance",
        "Average and RMS values of AC, phasor diagrams, power factor, active/reactive power, "
        "series RLC circuit and condition for resonance."
    )


def physics_lasers(request):
    return topic_placeholder(
        request,
        "physics",
        "Laser Physics",
        "Spontaneous and stimulated emission, population inversion, basic working of ruby/He–Ne laser, "
        "properties of laser beam and simple applications."
    )


def physics_optical_fiber(request):
    return topic_placeholder(
        request,
        "physics",
        "Optical Fiber & Waveguides",
        "Total internal reflection, numerical aperture, acceptance angle, propagation of light in optical fibers "
        "and use of fibers in communication and medicine."
    )


def physics_nuclear(request):
    return topic_placeholder(
        request,
        "physics",
        "Nuclear Physics & Radioactivity",
        "Nuclear composition, binding energy curve, radioactivity laws, half-life, nuclear reactions and basic ideas of fission & fusion."
    )


def physics_errors_data(request):
    return topic_placeholder(
        request,
        "physics",
        "Error & Data Analysis",
        "Types of errors (systematic, random), least count, absolute/relative error, propagation of errors and plotting graphs with best-fit line."
    )


def physics_instruments(request):
    return topic_placeholder(
        request,
        "physics",
        "Instruments & Measurements",
        "Working and least count of vernier caliper, screw gauge, spherometer, ammeter, voltmeter, galvanometer and meter bridge (conceptual level)."
    )

def physics_rotational(request):
    return topic_placeholder(
        request,
        "physics",
        "Rotational Motion & Rigid Bodies",
        "Torque, angular momentum, moment of inertia, rolling motion, rotational kinetic energy, parallel & perpendicular axis theorem."
    )


def physics_solids(request):
    return topic_placeholder(
        request,
        "physics",
        "Properties of Solids & Elasticity",
        "Stress-strain, Young’s modulus, bulk modulus, shear modulus, elastic potential energy, Hooke’s law and breaking stress."
    )


def physics_kinetic_theory(request):
    return topic_placeholder(
        request,
        "physics",
        "Kinetic Theory of Gases",
        "Assumptions of kinetic theory, pressure of an ideal gas, rms speed, mean free path, degrees of freedom and specific heats."
    )


def physics_semiconductors(request):
    return topic_placeholder(
        request,
        "physics",
        "Semiconductors & Electronics",
        "Intrinsic & extrinsic semiconductors, pn junction diode, Zener diode, transistor basics, rectifiers and simple logic gates."
    )


def physics_communication(request):
    return topic_placeholder(
        request,
        "physics",
        "Communication Systems",
        "Elements of communication system, analog & digital signals, modulation (AM, FM), bandwidth and basic idea of satellite communication."
    )


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


def physics_rotational(request):
    result = error = None
    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # Torque τ = I α
            if calc == "torque":
                I = float(request.POST["I"])
                a = float(request.POST["a"])
                torque = I * a
                result = f"Torque τ = {torque:.4f} N·m"

            # ω = ω₀ + αt
            elif calc == "omega":
                w0 = float(request.POST["w0"])
                a = float(request.POST["a"])
                t = float(request.POST["t"])
                w = w0 + a * t
                result = f"Angular velocity ω = {w:.4f} rad/s"

            # KE rotational: K = ½ I ω²
            elif calc == "ke_rot":
                I = float(request.POST["I"])
                w = float(request.POST["w"])
                K = 0.5 * I * w * w
                result = f"Rotational KE = {K:.4f} J"

        except Exception as e:
            error = f"Error: {e}"

    return render(request, "portal/physics_rotational.html",
                  {"active_subject": "physics", "result": result, "error": error})


def physics_solids(request):
    result = error = None
    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # σ = F/A
            if calc == "stress":
                F = float(request.POST["F"])
                A = float(request.POST["A"])
                stress = F / A
                result = f"Stress σ = {stress:.4f} N/m²"

            # ε = ΔL / L0
            elif calc == "strain":
                dL = float(request.POST["dL"])
                L0 = float(request.POST["L0"])
                strain = dL / L0
                result = f"Strain ε = {strain:.6f}"

            # Y = σ / ε
            elif calc == "young":
                sigma = float(request.POST["sigma"])
                eps = float(request.POST["eps"])
                Y = sigma / eps
                result = f"Young's modulus Y = {Y:.4e} N/m²"

        except Exception as e:
            error = f"Error: {e}"

    return render(request, "portal/physics_solids.html",
                  {"active_subject": "physics", "result": result, "error": error})


def physics_kinetic_theory(request):
    result = error = None
    R = 8.314
    if request.method == "POST":
        try:
            T = float(request.POST["T"])
            M = float(request.POST["M"])
            vrms = (3 * R * T / M) ** 0.5
            result = f"Root mean square speed vᵣₘₛ = {vrms:.4f} m/s"
        except Exception as e:
            error = f"Error: {e}"
    return render(request, "portal/physics_kinetic_theory.html",
                  {"active_subject": "physics", "result": result, "error": error})


def physics_semiconductors(request):
    result = error = None
    if request.method == "POST":
        try:
            I0 = float(request.POST["I0"])
            V = float(request.POST["V"])
            eta = float(request.POST["eta"])
            Vt = 0.026
            I = I0 * (math.exp(V / (eta * Vt)) - 1)
            result = f"Diode current I = {I:.4e} A"
        except Exception as e:
            error = f"Error: {e}"
    return render(request, "portal/physics_semiconductors.html",
                  {"active_subject": "physics", "result": result, "error": error})


def physics_communication(request):
    result = error = None
    if request.method == "POST":
        try:
            Amax = float(request.POST["Amax"])
            Amin = float(request.POST["Amin"])
            m = (Amax - Amin) / (Amax + Amin)
            result = f"Modulation Index m = {m:.4f}"
        except Exception as e:
            error = f"Error: {e}"
    return render(request, "portal/physics_communication.html",
                  {"active_subject": "physics", "result": result, "error": error})


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

# ===================== PHYSICS — NEW INTERACTIVE TOPICS =====================

def physics_units(request):
    """
    Units, Dimensions & Errors:
    - Simple length unit conversion
    - Time unit conversion
    - Percentage error from true & measured value
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # Length: m to km and cm
            if calc == "length":
                m = float(request.POST["m"])
                km = m / 1000.0
                cm = m * 100.0
                result = f"{m:.4f} m = {km:.6f} km = {cm:.2f} cm"

            # Time: seconds to minutes & hours
            elif calc == "time":
                s = float(request.POST["s"])
                minutes = s / 60.0
                hours = s / 3600.0
                result = f"{s:.4f} s = {minutes:.4f} min = {hours:.4f} h"

            # Percentage error
            elif calc == "per_error":
                true = float(request.POST["true"])
                meas = float(request.POST["meas"])
                if true == 0:
                    error = "True value cannot be zero for percentage error."
                else:
                    per = abs(meas - true) / abs(true) * 100
                    result = f"Percentage error = {per:.2f}%"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_units.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_ray_optics(request):
    """
    Ray Optics:
    - Lens/mirror formula: 1/f = 1/v - 1/u
    - Linear magnification: m = v/u
    - Power of lens: P = 1/f (in metres)
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            # Find v from u and f
            if calc == "lens_v":
                u = float(request.POST["u"])
                f = float(request.POST["f"])
                if (1.0/f + 1.0/u) == 0:
                    error = "Values give division by zero, change u or f."
                else:
                    v = 1.0 / (1.0/f + 1.0/u)
                    result = f"Image distance v = {v:.4f} m"

            # Magnification m = v/u
            elif calc == "magnification":
                u = float(request.POST["u"])
                v = float(request.POST["v"])
                if u == 0:
                    error = "Object distance u cannot be zero."
                else:
                    m = v / u
                    result = f"Magnification m = {m:.4f}"

            # Power P = 1/f (f in m, P in dioptre)
            elif calc == "power":
                f_cm = float(request.POST["f_cm"])
                if f_cm == 0:
                    error = "Focal length cannot be zero."
                else:
                    f_m = f_cm / 100.0
                    P = 1.0 / f_m
                    result = f"Power of lens P = {P:.2f} D"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_ray_optics.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_wave_optics(request):
    """
    Wave Optics:
    - YDSE fringe width: β = λ D / d
    - Path difference from fringe order: Δx = n β
    - Angular width of central maximum (single slit): θ = 2 λ / a
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "fringe_width":
                lam = float(request.POST["lam"])
                D = float(request.POST["D"])
                d = float(request.POST["d"])
                if d == 0:
                    error = "Slit separation d cannot be zero."
                else:
                    beta = lam * D / d
                    result = f"Fringe width β = {beta:.4e} m"

            elif calc == "path_diff":
                n = float(request.POST["n"])
                beta = float(request.POST["beta"])
                dx = n * beta
                result = f"Path difference Δx = {dx:.4e} m"

            elif calc == "single_slit":
                lam = float(request.POST["lam"])
                a = float(request.POST["a"])
                if a == 0:
                    error = "Slit width a cannot be zero."
                else:
                    theta = 2 * lam / a
                    result = f"Angular width θ ≈ {theta:.4e} rad"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_wave_optics.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_thermal_properties(request):
    """
    Thermal Properties of Matter:
    - Heat: Q = m s ΔT
    - Linear expansion: ΔL = α L ΔT
    - Latent heat: Q = m L
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "heat":
                m = float(request.POST["m"])
                s = float(request.POST["s"])
                dT = float(request.POST["dT"])
                Q = m * s * dT
                result = f"Heat Q = {Q:.4f} J"

            elif calc == "expansion":
                alpha = float(request.POST["alpha"])
                L = float(request.POST["L"])
                dT = float(request.POST["dT"])
                dL = alpha * L * dT
                result = f"Change in length ΔL = {dL:.6f} m"

            elif calc == "latent":
                m = float(request.POST["m"])
                Lh = float(request.POST["Lh"])
                Q = m * Lh
                result = f"Latent heat Q = {Q:.4f} J"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_thermal_properties.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_capacitance(request):
    """
    Electric Potential & Capacitance:
    - Q = C V
    - C (parallel plate) = ε A / d
    - Energy stored: U = ½ C V²
    """
    result = None
    error = None
    eps0 = 8.854e-12

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "q_cv":
                C = float(request.POST["C"])
                V = float(request.POST["V"])
                Q = C * V
                result = f"Charge Q = {Q:.4e} C"

            elif calc == "parallel_plate":
                rel = float(request.POST["rel_eps"])
                A = float(request.POST["A"])
                d = float(request.POST["d"])
                if d == 0:
                    error = "Plate separation d cannot be zero."
                else:
                    C = rel * eps0 * A / d
                    result = f"Capacitance C = {C:.4e} F"

            elif calc == "energy":
                C = float(request.POST["C"])
                V = float(request.POST["V"])
                U = 0.5 * C * V * V
                result = f"Energy stored U = {U:.4f} J"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_capacitance.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_magnetic_effects(request):
    """
    Magnetic Effects of Current:
    - B near long straight wire: B = μ₀ I / (2π r)
    - Force on conductor: F = B I L sinθ
    - Force on moving charge: F = q v B sinθ
    """
    result = None
    error = None
    mu0 = 4 * math.pi * 1e-7

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "wire_B":
                I = float(request.POST["I"])
                r = float(request.POST["r"])
                if r == 0:
                    error = "Distance r cannot be zero."
                else:
                    B = mu0 * I / (2 * math.pi * r)
                    result = f"Magnetic field B = {B:.4e} T"

            elif calc == "force_wire":
                B = float(request.POST["B"])
                I = float(request.POST["I"])
                L = float(request.POST["L"])
                theta = float(request.POST["theta"])
                F = B * I * L * math.sin(math.radians(theta))
                result = f"Force on conductor F = {F:.4e} N"

            elif calc == "force_charge":
                q = float(request.POST["q"])
                v = float(request.POST["v"])
                B = float(request.POST["B"])
                theta = float(request.POST["theta"])
                F = q * v * B * math.sin(math.radians(theta))
                result = f"Force on charge F = {F:.4e} N"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_magnetic_effects.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_emi_detailed(request):
    """
    Electromagnetic Induction (Detailed):
    - ε = - dΦ/dt
    - Rotating coil: ε_max = N B A ω
    - Inductor: energy U = ½ L I²
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "basic_emf":
                dphi = float(request.POST["dphi"])
                dt = float(request.POST["dt"])
                if dt == 0:
                    error = "Δt cannot be zero."
                else:
                    emf = - dphi / dt
                    result = f"Induced emf ε = {emf:.4e} V (sign from Lenz's law)"

            elif calc == "rotating_coil":
                N = float(request.POST["N"])
                B = float(request.POST["B"])
                A = float(request.POST["A"])
                w = float(request.POST["w"])
                emf = N * B * A * w
                result = f"Maximum emf ε_max = {emf:.4e} V"

            elif calc == "inductor_energy":
                L = float(request.POST["L"])
                I = float(request.POST["I"])
                U = 0.5 * L * I * I
                result = f"Energy stored in inductor U = {U:.4f} J"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_emi_detailed.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_ac_power(request):
    """
    AC Power & Resonance:
    - Reactances: X_L = 2π f L, X_C = 1/(2π f C)
    - Impedance: Z = √(R² + (X_L - X_C)²)
    - Resonant frequency: f₀ = 1 / (2π √(L C))
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "reactances":
                f = float(request.POST["f"])
                L = float(request.POST["L"])
                C = float(request.POST["C"])
                if f <= 0:
                    error = "Frequency must be > 0."
                else:
                    XL = 2 * math.pi * f * L
                    XC = 0.0 if C <= 0 else 1.0 / (2 * math.pi * f * C)
                    result = f"X_L = {XL:.4f} Ω,  X_C = {XC:.4f} Ω"

            elif calc == "impedance":
                R = float(request.POST["R"])
                XL = float(request.POST["XL"])
                XC = float(request.POST["XC"])
                Z = math.sqrt(R * R + (XL - XC) ** 2)
                result = f"Impedance Z = {Z:.4f} Ω"

            elif calc == "resonance":
                L = float(request.POST["L"])
                C = float(request.POST["C"])
                if L <= 0 or C <= 0:
                    error = "L and C must both be > 0."
                else:
                    f0 = 1.0 / (2 * math.pi * math.sqrt(L * C))
                    result = f"Resonant frequency f₀ = {f0:.4f} Hz"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_ac_power.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_lasers(request):
    """
    Laser Physics:
    - Photon energy: E = h c / λ
    - Number of photons from power & time: N = (P t) / E
    """
    result = None
    error = None
    h = 6.626e-34
    c = 3e8

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "photon_energy":
                lam = float(request.POST["lam"])
                if lam <= 0:
                    error = "Wavelength λ must be > 0."
                else:
                    E = h * c / lam
                    result = f"Photon energy E = {E:.4e} J"

            elif calc == "num_photons":
                lam = float(request.POST["lam"])
                P = float(request.POST["P"])
                t = float(request.POST["t"])
                if lam <= 0:
                    error = "Wavelength λ must be > 0."
                else:
                    E = h * c / lam
                    N = P * t / E
                    result = f"Approx. number of photons N ≈ {N:.4e}"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_lasers.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_optical_fiber(request):
    """
    Optical Fiber & Waveguides:
    - Numerical aperture: NA = √(n₁² - n₂²)
    - Acceptance angle: θₐ = sin⁻¹(NA)
    - Critical angle: sin c = n₂ / n₁
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "na":
                n1 = float(request.POST["n1"])
                n2 = float(request.POST["n2"])
                if n1 <= n2:
                    error = "Core index n₁ must be > cladding index n₂."
                else:
                    NA = math.sqrt(n1 * n1 - n2 * n2)
                    result = f"Numerical aperture NA = {NA:.4f}"

            elif calc == "accept_angle":
                NA = float(request.POST["NA"])
                if NA < 0 or NA > 1:
                    error = "NA should be between 0 and 1."
                else:
                    theta = math.degrees(math.asin(NA))
                    result = f"Acceptance angle θₐ ≈ {theta:.2f}°"

            elif calc == "critical_angle":
                n1 = float(request.POST["n1"])
                n2 = float(request.POST["n2"])
                if n1 == 0 or n2 > n1:
                    error = "Need n₂ ≤ n₁ and n₁ > 0."
                else:
                    cang = math.degrees(math.asin(n2 / n1))
                    result = f"Critical angle c ≈ {cang:.2f}°"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_optical_fiber.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_nuclear(request):
    """
    Nuclear Physics & Radioactivity:
    - Decay law: N = N₀ e^{-λ t}
    - Activity A = λ N
    - Relation of half-life: λ = ln 2 / T₁/₂
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "decay":
                N0 = float(request.POST["N0"])
                lam = float(request.POST["lam"])
                t = float(request.POST["t"])
                N = N0 * math.exp(-lam * t)
                result = f"N(t) = {N:.4e}"

            elif calc == "activity":
                lam = float(request.POST["lam"])
                N = float(request.POST["N"])
                A = lam * N
                result = f"Activity A = {A:.4e} decays/s"

            elif calc == "lambda_from_t12":
                t12 = float(request.POST["t12"])
                if t12 <= 0:
                    error = "Half-life must be > 0."
                else:
                    lam = math.log(2) / t12
                    result = f"Decay constant λ = {lam:.4e} s⁻¹"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_nuclear.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_errors_data(request):
    """
    Error & Data Analysis:
    - Mean of measurements
    - Absolute & percentage error of one reading w.r.t true
    - Standard deviation of measurements
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "mean":
                raw = request.POST["values"]
                nums = [float(x) for x in raw.replace(" ", "").split(",") if x != ""]
                if not nums:
                    error = "Enter at least one value."
                else:
                    mean = sum(nums) / len(nums)
                    result = f"Mean of measurements = {mean:.4f}"

            elif calc == "single_error":
                true = float(request.POST["true"])
                meas = float(request.POST["meas"])
                if true == 0:
                    error = "True value cannot be zero."
                else:
                    abs_err = abs(meas - true)
                    per = abs_err / abs(true) * 100
                    result = f"Absolute error = {abs_err:.4f}, Percentage error = {per:.2f}%"

            elif calc == "stddev":
                raw = request.POST["values"]
                nums = [float(x) for x in raw.replace(" ", "").split(",") if x != ""]
                if not nums:
                    error = "Enter at least one value."
                else:
                    n = len(nums)
                    mean = sum(nums) / n
                    var = sum((x - mean) ** 2 for x in nums) / n
                    sd = math.sqrt(var)
                    result = f"Standard deviation σ = {sd:.4f}"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_errors_data.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


def physics_instruments(request):
    """
    Instruments & Measurements:
    - Least count of vernier: LC = 1 MSD - 1 VSD
    - Least count of screw gauge: LC = pitch / no. of divisions
    - Percentage error from least count and measured magnitude
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "vernier_lc":
                msd = float(request.POST["msd"])
                vsd = float(request.POST["vsd"])
                if vsd == 0:
                    error = "Vernier division value cannot be zero."
                else:
                    lc = abs(msd - vsd)
                    result = f"Approx. vernier least count ≈ {lc:.4f} (same unit as scale)"

            elif calc == "screw_lc":
                pitch = float(request.POST["pitch"])
                divs = float(request.POST["divs"])
                if divs == 0:
                    error = "Number of divisions cannot be zero."
                else:
                    lc = pitch / divs
                    result = f"Screw gauge least count = {lc:.6f} (same unit as pitch)"

            elif calc == "lc_error":
                lc = float(request.POST["lc"])
                value = float(request.POST["value"])
                if value == 0:
                    error = "Measured value cannot be zero."
                else:
                    per = lc / value * 100
                    result = f"Approx. percentage error due to instrument = {per:.2f}%"

        except (KeyError, ValueError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/physics_instruments.html",
        {"active_subject": "physics", "result": result, "error": error},
    )


# ===================== MATHS — TOPIC PAGES + CALCULATORS =====================

# ---------- 1. Syllabus-style topic pages (mostly UI-only shells) ----------

def maths_real_numbers(request):
    """
    Real Numbers:
    Number system, Euclid’s division lemma, HCF/LCM etc.
    (Main content is in maths_real_numbers.html)
    """
    return render(
        request,
        "portal/maths_real_numbers.html",
        {"active_subject": "maths"},
    )


def maths_polynomials(request):
    """
    Polynomials:
    Yeh tumhara detailed 'Maths – Part 1' page hai with multiple calculators.
    """
    return render(
        request,
        "portal/maths_polynomials.html",
        {"active_subject": "maths"},
    )


def maths_linear_equations(request):
    """
    Linear Equations in Two Variables:
    Graph / intersection etc. – main UI in template.
    """
    return render(
        request,
        "portal/maths_linear_equations.html",
        {"active_subject": "maths"},
    )


def maths_quadratic_equations(request):
    """
    Quadratic Equations topic page.
    """
    return render(
        request,
        "portal/maths_quadratic_equations.html",
        {"active_subject": "maths"},
    )


# NOTE: maths_ap को हम नीचे interactive calculator की तरह ही use करेंगे,
# लेकिन url name वही रहेगा: 'maths_ap'.


# ---------- 2. Coordinate Geometry (Basics) ----------

def maths_coordinate_basics(request):
    """
    Coordinate Geometry (Basics):
    - Distance between two points
    - Internal division (section formula) of a line segment
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            x1 = float(request.POST.get("x1", 0))
            y1 = float(request.POST.get("y1", 0))
            x2 = float(request.POST.get("x2", 0))
            y2 = float(request.POST.get("y2", 0))

            if calc == "distance":
                d = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
                result = f"Distance between points = {d:.4f}"

            elif calc == "section":
                m = float(request.POST.get("m", 0))
                n = float(request.POST.get("n", 0))
                if m + n == 0:
                    error = "m + n cannot be zero."
                else:
                    X = (m * x2 + n * x1) / (m + n)
                    Y = (m * y2 + n * y1) / (m + n)
                    result = f"Internal division point = ({X:.4f}, {Y:.4f})"
            else:
                error = "Unknown calculation type."
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_coordinate_basics.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


# ---------- 3. Trigonometry (Basics + Applications) ----------

def maths_trig_basic(request):
    """
    Trigonometry – Basics:
    - sin, cos, tan of an angle in degrees
    - Degree ↔ Radian conversion
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "values":
                angle_deg = float(request.POST.get("angle", 0))
                ang_rad = math.radians(angle_deg)
                s = math.sin(ang_rad)
                c = math.cos(ang_rad)
                if abs(c) < 1e-12:
                    t = "undefined"
                else:
                    t = f"{math.tan(ang_rad):.4f}"
                result = (
                    f"sin({angle_deg}°) = {s:.4f}, "
                    f"cos({angle_deg}°) = {c:.4f}, "
                    f"tan({angle_deg}°) = {t}"
                )

            elif calc == "deg2rad":
                deg = float(request.POST.get("deg", 0))
                rad = math.radians(deg)
                result = f"{deg}° = {rad:.6f} rad"

            elif calc == "rad2deg":
                rad = float(request.POST.get("rad", 0))
                deg = math.degrees(rad)
                result = f"{rad} rad = {deg:.4f}°"

            else:
                error = "Unknown calculation type."
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_trig_basic.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_trig_applications(request):
    """
    Trigonometry – Applications:
    - Height from distance and angle (tan θ = h/d)
    - Distance from height and angle
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            angle_deg = float(request.POST.get("angle", 0))
            ang_rad = math.radians(angle_deg)

            if calc == "height_from_distance":
                d = float(request.POST.get("d", 0))
                h = d * math.tan(ang_rad)
                result = f"Height = {h:.4f} units"

            elif calc == "distance_from_height":
                h = float(request.POST.get("h", 0))
                tan_val = math.tan(ang_rad)
                if abs(tan_val) < 1e-12:
                    error = "Angle too small; tan(θ) ≈ 0 → distance becomes infinite."
                else:
                    d = h / tan_val
                    result = f"Horizontal distance = {d:.4f} units"
            else:
                error = "Unknown calculation type."

        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_trig_applications.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


# ---------- 4. Circles & Areas Related to Circles ----------

def maths_circles(request):
    """
    Circles & Tangents:
    - Circumference and area of a circle
    - Length of tangent from an external point: √(d² − r²)
    """
    result = None
    error = None
    pi = math.pi

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "circle":
                r = float(request.POST.get("r", 0))
                if r < 0:
                    error = "Radius cannot be negative."
                else:
                    C = 2 * pi * r
                    A = pi * r * r
                    result = (
                        f"Circumference = {C:.4f} units, "
                        f"Area = {A:.4f} square units"
                    )

            elif calc == "tangent":
                d = float(request.POST.get("d", 0))
                r = float(request.POST.get("r", 0))
                if d <= r:
                    error = "External point must be outside the circle (d > r)."
                else:
                    L = math.sqrt(d * d - r * r)
                    result = f"Length of tangent from external point = {L:.4f} units"

            else:
                error = "Unknown calculation type."
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_circles.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_area_circles(request):
    """
    Areas Related to Circles:
    - Area of sector: (θ/360) π r²
    - Arc length: (θ/360) 2 π r
    """
    result = None
    error = None
    pi = math.pi

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            r = float(request.POST.get("r", 0))
            theta = float(request.POST.get("theta", 0))  # in degrees

            if r < 0:
                error = "Radius cannot be negative."
            else:
                if calc == "sector_area":
                    A = (theta / 360.0) * pi * r * r
                    result = f"Area of sector = {A:.4f} square units"

                elif calc == "arc_length":
                    L = (theta / 360.0) * 2 * pi * r
                    result = f"Arc length = {L:.4f} units"
                else:
                    error = "Unknown calculation type."
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_area_circles.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


# ---------- 5. Surface Area & Volume ----------

def maths_surface_volume(request):
    """
    Surface Area & Volume:
    - TSA & volume of cuboid (l, b, h)
    - CSA/TSA & volume of cylinder (r, h)
    - TSA & volume of cone (r, h)
    """
    result = None
    error = None
    pi = math.pi

    if request.method == "POST":
        shape = request.POST.get("shape")

        try:
            # CUBOID
            if shape == "cuboid":
                l = float(request.POST.get("l", 0))
                b = float(request.POST.get("b", 0))
                h = float(request.POST.get("h", 0))
                if l < 0 or b < 0 or h < 0:
                    error = "Dimensions cannot be negative."
                else:
                    tsa = 2 * (l * b + b * h + h * l)
                    vol = l * b * h
                    result = (
                        f"Cuboid TSA = {tsa:.4f} square units, "
                        f"Volume = {vol:.4f} cubic units"
                    )

            # CYLINDER
            elif shape == "cylinder":
                r = float(request.POST.get("r", 0))
                h = float(request.POST.get("h", 0))
                if r < 0 or h < 0:
                    error = "Dimensions cannot be negative."
                else:
                    csa = 2 * pi * r * h
                    tsa = csa + 2 * pi * r * r
                    vol = pi * r * r * h
                    result = (
                        f"Cylinder CSA = {csa:.4f}, TSA = {tsa:.4f} square units, "
                        f"Volume = {vol:.4f} cubic units"
                    )

            # CONE
            elif shape == "cone":
                r = float(request.POST.get("r", 0))
                h = float(request.POST.get("h", 0))
                if r < 0 or h < 0:
                    error = "Dimensions cannot be negative."
                else:
                    l = math.sqrt(r * r + h * h)  # slant height
                    csa = pi * r * l
                    tsa = csa + pi * r * r
                    vol = (1.0 / 3.0) * pi * r * r * h
                    result = (
                        f"Cone CSA = {csa:.4f}, TSA = {tsa:.4f} square units, "
                        f"Volume = {vol:.4f} cubic units"
                    )

            else:
                error = "Unknown shape selected."

        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_surface_volume.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


# ---------- 6. Relations & Functions, Matrices, Determinants ----------

def maths_relations_functions(request):
    """
    Relations & Functions:
    - Evaluate linear function f(x) = a x + b
    - Inverse of linear function (if a ≠ 0)
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            a = float(request.POST.get("a", 0))
            b = float(request.POST.get("b", 0))

            if calc == "value":
                x = float(request.POST.get("x", 0))
                fx = a * x + b
                result = f"f({x}) = {fx:.4f}"

            elif calc == "inverse":
                if a == 0:
                    error = "Inverse does not exist if a = 0 (not one–one)."
                else:
                    result = (
                        "Inverse function: f⁻¹(y) = (y − b) / a\n"
                        f"Here: f⁻¹(y) = (y − {b}) / {a}"
                    )
            else:
                error = "Unknown calculation type."
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_relations_functions.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_matrices(request):
    """
    Matrices:
    - Addition of two 2×2 matrices
    - Scalar multiplication of a 2×2 matrix
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")

        try:
            if calc == "add_2x2":
                a11 = float(request.POST.get("a11", 0))
                a12 = float(request.POST.get("a12", 0))
                a21 = float(request.POST.get("a21", 0))
                a22 = float(request.POST.get("a22", 0))

                b11 = float(request.POST.get("b11", 0))
                b12 = float(request.POST.get("b12", 0))
                b21 = float(request.POST.get("b21", 0))
                b22 = float(request.POST.get("b22", 0))

                c11 = a11 + b11
                c12 = a12 + b12
                c21 = a21 + b21
                c22 = a22 + b22

                result = (
                    f"Result C = A + B = [[{c11:.2f}, {c12:.2f}], "
                    f"[{c21:.2f}, {c22:.2f}]]"
                )

            elif calc == "scalar_2x2":
                k = float(request.POST.get("k", 0))
                a11 = float(request.POST.get("a11", 0))
                a12 = float(request.POST.get("a12", 0))
                a21 = float(request.POST.get("a21", 0))
                a22 = float(request.POST.get("a22", 0))

                c11 = k * a11
                c12 = k * a12
                c21 = k * a21
                c22 = k * a22

                result = (
                    f"kA = [[{c11:.2f}, {c12:.2f}], "
                    f"[{c21:.2f}, {c22:.2f}]]"
                )
            else:
                error = "Unknown calculation type."
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_matrices.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_determinants(request):
    """
    Determinants:
    - 2×2 determinant: |a b; c d|
    - 3×3 determinant using Sarrus rule
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")

        try:
            if calc == "det2":
                a = float(request.POST.get("a", 0))
                b = float(request.POST.get("b", 0))
                c = float(request.POST.get("c", 0))
                d = float(request.POST.get("d", 0))
                det = a * d - b * c
                result = f"Determinant = {det:.4f}"

            elif calc == "det3":
                a = float(request.POST.get("a", 0))
                b = float(request.POST.get("b", 0))
                c = float(request.POST.get("c", 0))
                d = float(request.POST.get("d", 0))
                e = float(request.POST.get("e", 0))
                f = float(request.POST.get("f", 0))
                g = float(request.POST.get("g", 0))
                h = float(request.POST.get("h", 0))
                i = float(request.POST.get("i", 0))

                det = (
                    a * (e * i - f * h)
                    - b * (d * i - f * g)
                    + c * (d * h - e * g)
                )
                result = f"Determinant = {det:.4f}"
            else:
                error = "Unknown calculation type."
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_determinants.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


# ---------- 7. Limits, Differentiation, Applications, Integration ----------

def maths_limits_continuity(request):
    """
    Limits & Continuity:
    - Approximate left & right limits of f(x) = ax² + bx + c at x = x₀ using step h
    """
    result = None
    error = None

    if request.method == "POST":
        try:
            a = float(request.POST.get("a", 0))
            b = float(request.POST.get("b", 0))
            c = float(request.POST.get("c", 0))
            x0 = float(request.POST.get("x0", 0))
            h = float(request.POST.get("h", 0))

            if h <= 0:
                error = "h must be positive (small step)."
            else:
                def f(x):
                    return a * x * x + b * x + c

                left = f(x0 - h)
                right = f(x0 + h)
                at = f(x0)

                result = (
                    f"f(x) = {a}x² + {b}x + {c}\n"
                    f"Left approx f(x₀ − h) = {left:.4f}\n"
                    f"Right approx f(x₀ + h) = {right:.4f}\n"
                    f"f(x₀) = {at:.4f}"
                )
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_limits_continuity.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_differentiation(request):
    """
    Differentiation:
    - f(x) = ax² + bx + c → f'(x) = 2ax + b
    """
    result = None
    error = None

    if request.method == "POST":
        try:
            a = float(request.POST.get("a", 0))
            b = float(request.POST.get("b", 0))
            c = float(request.POST.get("c", 0))  # not used
            x = float(request.POST.get("x", 0))

            fp = 2 * a * x + b
            result = f"f'(x) at x = {x} is {fp:.4f}"
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_differentiation.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_app_derivatives(request):
    """
    Applications of Derivatives:
    - Slope of tangent for y = ax² + bx + c at x₀
    - Nature of extremum (max/min) for quadratic
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            a = float(request.POST.get("a", 0))
            b = float(request.POST.get("b", 0))
            c = float(request.POST.get("c", 0))

            if calc == "slope":
                x0 = float(request.POST.get("x0", 0))
                slope = 2 * a * x0 + b
                result = f"Slope of tangent at x = {x0} is {slope:.4f}"

            elif calc == "extremum":
                if a == 0:
                    error = "Not a quadratic (a = 0), cannot decide max/min."
                else:
                    xv = -b / (2 * a)
                    yv = a * xv * xv + b * xv + c
                    nature = "Minimum" if a > 0 else "Maximum"
                    result = f"{nature} value at x = {xv:.4f} is y = {yv:.4f}"
            else:
                error = "Unknown calculation type."
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_app_derivatives.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_integration(request):
    """
    Integration:
    - ∫ (ax² + bx + c) dx = (a/3)x³ + (b/2)x² + cx + C
    - Evaluate integral at a particular x
    """
    result = None
    error = None

    if request.method == "POST":
        try:
            a = float(request.POST.get("a", 0))
            b = float(request.POST.get("b", 0))
            c = float(request.POST.get("c", 0))
            x = float(request.POST.get("x", 0))

            val = (a / 3.0) * x ** 3 + (b / 2.0) * x ** 2 + c * x
            result = (
                "∫(ax² + bx + c) dx = (a/3)x³ + (b/2)x² + cx + C\n"
                f"At x = {x}, integral value (without C) = {val:.4f}"
            )
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_integration.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_app_integrals(request):
    """
    Applications of Integrals:
    - Area under line y = mx + c from x = x1 to x = x2
    """
    result = None
    error = None

    if request.method == "POST":
        try:
            m = float(request.POST.get("m", 0))
            c0 = float(request.POST.get("c0", 0))
            x1 = float(request.POST.get("x1", 0))
            x2 = float(request.POST.get("x2", 0))

            if x1 == x2:
                error = "x1 and x2 must be different."
            else:
                def F(x):
                    return (m / 2.0) * x ** 2 + c0 * x

                area = F(x2) - F(x1)
                if area < 0:
                    area = -area
                result = f"Area under y = {m}x + {c0} from {x1} to {x2} ≈ {area:.4f}"
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_app_integrals.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_diff_eq(request):
    """
    Differential Equations:
    - Solve dy/dx = k y with y(0) = y0 → y(x) = y0 e^{kx}
    """
    result = None
    error = None

    if request.method == "POST":
        try:
            k = float(request.POST.get("k", 0))
            y0 = float(request.POST.get("y0", 0))
            x = float(request.POST.get("x", 0))

            y = y0 * math.exp(k * x)
            result = f"Solution: y(x) = {y0} e^({k}x), so y({x}) = {y:.4f}"
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_diff_eq.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


# ---------- 8. Vector Algebra, 3D Geometry, Linear Programming ----------

def maths_vectors(request):
    """
    Vector Algebra:
    - Magnitude of 3D vector
    - Dot product and cross product of two 3D vectors
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")

        try:
            ax = float(request.POST.get("ax", 0))
            ay = float(request.POST.get("ay", 0))
            az = float(request.POST.get("az", 0))

            if calc == "mag":
                mag = math.sqrt(ax * ax + ay * ay + az * az)
                result = f"|a| = {mag:.4f}"

            else:
                bx = float(request.POST.get("bx", 0))
                by = float(request.POST.get("by", 0))
                bz = float(request.POST.get("bz", 0))

                if calc == "dot":
                    dot = ax * bx + ay * by + az * bz
                    result = f"a · b = {dot:.4f}"

                elif calc == "cross":
                    cx = ay * bz - az * by
                    cy = az * bx - ax * bz
                    cz = ax * by - ay * bx
                    result = f"a × b = ({cx:.4f}, {cy:.4f}, {cz:.4f})"
                else:
                    error = "Unknown calculation type."
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_vectors.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_3d_geometry(request):
    """
    3D Geometry:
    - Distance between two points in 3D
    - Distance of a point from origin
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            x1 = float(request.POST.get("x1", 0))
            y1 = float(request.POST.get("y1", 0))
            z1 = float(request.POST.get("z1", 0))

            if calc == "origin":
                d = math.sqrt(x1 * x1 + y1 * y1 + z1 * z1)
                result = f"Distance from origin = {d:.4f}"

            elif calc == "points":
                x2 = float(request.POST.get("x2", 0))
                y2 = float(request.POST.get("y2", 0))
                z2 = float(request.POST.get("z2", 0))
                d = math.sqrt(
                    (x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2
                )
                result = f"Distance between points = {d:.4f}"
            else:
                error = "Unknown calculation type."

        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_3d_geometry.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


def maths_linear_programming(request):
    """
    Linear Programming (very basic helper):
    - For two lines a1x + b1y = c1 and a2x + b2y = c2,
      find intersection and evaluate Z = px + qy there.
    """
    result = None
    error = None

    if request.method == "POST":
        try:
            a1 = float(request.POST.get("a1", 0))
            b1 = float(request.POST.get("b1", 0))
            c1 = float(request.POST.get("c1", 0))

            a2 = float(request.POST.get("a2", 0))
            b2 = float(request.POST.get("b2", 0))
            c2 = float(request.POST.get("c2", 0))

            p = float(request.POST.get("p", 0))
            q = float(request.POST.get("q", 0))

            D = a1 * b2 - a2 * b1
            if D == 0:
                error = "Lines are parallel or coincident; no unique intersection."
            else:
                Dx = c1 * b2 - c2 * b1
                Dy = a1 * c2 - a2 * c1
                x = Dx / D
                y = Dy / D
                Z = p * x + q * y
                result = (
                    f"Intersection at (x, y) = ({x:.4f}, {y:.4f})\n"
                    f"Objective Z = {p}x + {q}y = {Z:.4f}"
                )
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_linear_programming.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


# ---------- 9. Simple calculator tools: AP, Arithmetic, Algebra, Geometry, etc. ----------

def maths_ap(request):
    """
    Arithmetic Progressions (AP):
    - n-th term: aₙ = a + (n−1)d
    - Sum of n terms: Sₙ = n/2 [2a + (n−1)d]
    """
    result = None
    error = None

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            a = float(request.POST.get("a", 0))
            d = float(request.POST.get("d", 0))
            n = float(request.POST.get("n", 0))

            if n <= 0:
                error = "n must be positive."
            else:
                if calc == "nth":
                    an = a + (n - 1) * d
                    result = f"n-th term aₙ = {an:.4f}"
                elif calc == "sum":
                    Sn = (n / 2.0) * (2 * a + (n - 1) * d)
                    result = f"Sum of first {int(n)} terms Sₙ = {Sn:.4f}"
                else:
                    error = "Unknown calculation type."
        except (ValueError, TypeError) as e:
            error = f"Input error: {e}"

    return render(
        request,
        "portal/maths_ap.html",
        {"active_subject": "maths", "result": result, "error": error},
    )


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
    Coordinate Geometry (tool):
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


## ===================== CHEMISTRY — Full syllabus placeholder pages =====================

def chem_stoichiometry(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Stoichiometry & Mole Concept",
        "Mole concept, molar mass, percentage composition, limiting reagent and basic stoichiometric calculations."
    )


def chem_states_matter(request):
    return topic_placeholder(
        request,
        "chemistry",
        "States of Matter (Gas, Liquid & Solid)",
        "Gas laws, ideal gas equation, kinetic theory, intermolecular forces and types of solids."
    )


def chem_periodic_table(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Periodic Table & Periodicity",
        "Modern periodic table, classification, periodic trends like atomic radius, ionization enthalpy and electronegativity."
    )


def chem_equilibrium(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Chemical & Ionic Equilibrium",
        "Law of mass action, equilibrium constant, Le Chatelier’s principle, ionic equilibria, pH and solubility product."
    )


def chem_redox(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Redox Reactions",
        "Concept of oxidation and reduction, oxidation number, balancing redox reactions and applications."
    )


def chem_electrochemistry(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Electrochemistry",
        "Galvanic cells, EMF, Nernst equation, conductance, electrolytic cells and Faraday’s laws of electrolysis."
    )


def chem_chemical_kinetics(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Chemical Kinetics",
        "Rate of reaction, rate law, order and molecularity, integrated rate equations and Arrhenius equation."
    )


def chem_surface_chemistry(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Surface Chemistry",
        "Adsorption, catalysts, colloids, emulsions and their applications in daily life and industry."
    )


def chem_hydrogen(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Hydrogen & its Compounds",
        "Position of hydrogen, isotopes, hydrides, water, heavy water and hydrogen peroxide – preparation and properties."
    )


def chem_s_block(request):
    return topic_placeholder(
        request,
        "chemistry",
        "s-Block Elements",
        "General characteristics of alkali and alkaline earth metals, important compounds and biological significance."
    )


def chem_p_block(request):
    return topic_placeholder(
        request,
        "chemistry",
        "p-Block Elements",
        "Group 13–18 elements, important oxides, halides, oxyacids and their everyday applications."
    )


def chem_d_f_block(request):
    return topic_placeholder(
        request,
        "chemistry",
        "d- & f-Block Elements",
        "Transition and inner-transition elements, coloured ions, variable oxidation states and complex formation."
    )


def chem_coordination(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Coordination Compounds",
        "Ligands, coordination number, IUPAC nomenclature, isomerism and basic bonding ideas in coordination compounds."
    )


def chem_metallurgy(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Metallurgy & Qualitative Analysis",
        "General principles of extraction of metals, concentration, roasting, reduction and basic qualitative analysis."
    )


def chem_environmental(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Environmental Chemistry",
        "Types of pollution, greenhouse effect, ozone depletion, water and soil pollution and green chemistry."
    )


def chem_organic_basics(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Basic Concepts of Organic & IUPAC",
        "Classification of organic compounds, functional groups, IUPAC nomenclature, isomerism and electronic effects."
    )


def chem_hydrocarbons(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Hydrocarbons (Alkane, Alkene, Alkyne)",
        "Preparation, properties and reactions of alkanes, alkenes, alkynes and aromatic hydrocarbons."
    )


def chem_haloalkanes(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Haloalkanes & Haloarenes",
        "Nature of C–X bond, classification, nucleophilic substitution and elimination reactions and environmental impact."
    )


def chem_alcohols(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Alcohols, Phenols & Ethers",
        "Preparation, properties and important reactions of alcohols, phenols and ethers."
    )


def chem_aldehydes(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Aldehydes, Ketones & Carboxylic Acids",
        "Nomenclature, structure, nucleophilic addition reactions, acidity and important name reactions."
    )


def chem_amines(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Amines",
        "Classification, basic character, methods of preparation, reactions and diazonium salts."
    )


def chem_biomolecules(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Biomolecules",
        "Carbohydrates, proteins, vitamins and nucleic acids – basic structure and biological roles."
    )


def chem_polymers(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Polymers",
        "Classification, addition and condensation polymers, examples like PVC, teflon, bakelite and biodegradable polymers."
    )


def chem_everyday_life(request):
    return topic_placeholder(
        request,
        "chemistry",
        "Chemistry in Everyday Life",
        "Drugs, medicines, food additives, cleansing agents and simple chemical principles used in daily life."
    )


# ---------- CHEMISTRY: helper for syllabus topic pages ----------

def _chem_topic_page(request, title, description):
    """
    Reuse the common topic_placeholder template but with Chemistry subject.
    This keeps design exactly same as Physics/Maths placeholder pages.
    """
    return topic_placeholder(
        request,
        subject="chemistry",
        topic=title,
        description=description,
    )

def chem_stoichiometry(request):
    return _chem_topic_page(
        request,
        "Stoichiometry & Mole Concept",
        "Mole concept, Avogadro number, empirical and molecular formula, "
        "stoichiometric calculations with mass, volume and moles, limiting reagent and percentage yield."
    )


def chem_states_matter(request):
    return _chem_topic_page(
        request,
        "States of Matter (Gas, Liquid & Solid)",
        "Ideal and real gases, gas laws, kinetic theory basics, intermolecular forces, "
        "liquids and solids – properties and simple numerical ideas."
    )


def chem_periodic_table(request):
    return _chem_topic_page(
        request,
        "Periodic Table & Periodicity",
        "Modern periodic table, classification of elements, trends in atomic/ionic size, "
        "ionisation enthalpy, electron gain enthalpy, electronegativity and valency."
    )


def chem_equilibrium(request):
    return _chem_topic_page(
        request,
        "Chemical & Ionic Equilibrium",
        "Dynamic equilibrium, law of mass action, equilibrium constant, Le Chatelier principle, "
        "pH, buffer solutions, solubility product and common ion effect."
    )


def chem_redox(request):
    return _chem_topic_page(
        request,
        "Redox Reactions",
        "Oxidation–reduction concepts, oxidation number, balancing redox equations, "
        "redox titrations and applications in chemistry."
    )


def chem_electrochemistry(request):
    return _chem_topic_page(
        request,
        "Electrochemistry",
        "Electrolytic and galvanic cells, EMF of cell, Nernst equation, standard electrode potential, "
        "electrochemical series and simple batteries."
    )


def chem_chemical_kinetics(request):
    return _chem_topic_page(
        request,
        "Chemical Kinetics",
        "Rate of reaction, rate law and order, integrated rate expressions for simple reactions, "
        "Arrhenius equation and factors affecting rate."
    )


def chem_surface_chemistry(request):
    return _chem_topic_page(
        request,
        "Surface Chemistry",
        "Adsorption, catalysis, colloids, emulsions and applications of surface phenomena in daily life."
    )

def chem_hydrogen(request):
    return _chem_topic_page(
        request,
        "Hydrogen & its Compounds",
        "Position of hydrogen, isotopes, physical and chemical properties, "
        "hydrides, water, hydrogen peroxide and their uses."
    )


def chem_s_block(request):
    return _chem_topic_page(
        request,
        "s-Block Elements",
        "Group 1 and 2 elements – electronic configuration, trends, important compounds "
        "like NaCl, NaOH, CaO, CaCO₃ and their applications."
    )


def chem_p_block(request):
    return _chem_topic_page(
        request,
        "p-Block Elements",
        "Group 13–18 overview, important trends, and main compounds of boron, carbon, nitrogen, "
        "oxygen, halogens and noble gases."
    )


def chem_d_f_block(request):
    return _chem_topic_page(
        request,
        "d- & f-Block Elements",
        "Transition and inner-transition elements, variable oxidation states, colors, "
        "magnetic properties and important alloys."
    )


def chem_coordination(request):
    return _chem_topic_page(
        request,
        "Coordination Compounds",
        "Ligands, coordination number, Werner theory, nomenclature, isomerism and basic applications "
        "like cisplatin, hemoglobin, etc."
    )


def chem_metallurgy(request):
    return _chem_topic_page(
        request,
        "Metallurgy & Qualitative Analysis",
        "General principles of extraction of metals, concentration of ores, "
        "reduction, refining and a short overview of qualitative salt analysis."
    )


def chem_environmental(request):
    return _chem_topic_page(
        request,
        "Environmental Chemistry",
        "Pollution of air, water and soil, ozone depletion, greenhouse effect, "
        "acid rain and strategies for control of pollution."
    )
def chem_organic_basics(request):
    return _chem_topic_page(
        request,
        "Basic Concepts of Organic & IUPAC",
        "Classification of organic compounds, types of bonding and isomerism, "
        "IUPAC nomenclature and basics of reaction mechanisms."
    )


def chem_hydrocarbons(request):
    return _chem_topic_page(
        request,
        "Hydrocarbons (Alkane, Alkene, Alkyne)",
        "Preparation, properties and reactions of alkanes, alkenes and alkynes, "
        "including important named reactions and uses."
    )


def chem_haloalkanes(request):
    return _chem_topic_page(
        request,
        "Haloalkanes & Haloarenes",
        "Nomenclature, nature of C–X bond, substitution and elimination reactions, "
        "and environmental effects of halo-compounds."
    )


def chem_alcohols(request):
    return _chem_topic_page(
        request,
        "Alcohols, Phenols & Ethers",
        "Preparation and properties of alcohols, phenols and ethers with focus on "
        "important reactions and everyday applications."
    )


def chem_aldehydes(request):
    return _chem_topic_page(
        request,
        "Aldehydes, Ketones & Carboxylic Acids",
        "Structure, preparation and characteristic reactions including nucleophilic addition "
        "and acidity of carboxylic acids."
    )


def chem_amines(request):
    return _chem_topic_page(
        request,
        "Amines",
        "Classification, basicity, preparation and important reactions of aliphatic and aromatic amines."
    )


def chem_biomolecules(request):
    return _chem_topic_page(
        request,
        "Biomolecules",
        "Carbohydrates, proteins, vitamins and nucleic acids – basic structures and biological roles."
    )


def chem_polymers(request):
    return _chem_topic_page(
        request,
        "Polymers",
        "Types of polymers (addition, condensation, natural, synthetic), "
        "common examples and their uses in daily life."
    )


def chem_everyday_life(request):
    return _chem_topic_page(
        request,
        "Chemistry in Everyday Life",
        "Drugs, soaps and detergents, food additives, and other chemical substances used in daily life."
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







# ===================== CHEMISTRY — Interactive calculator pages =====================

def chem_atomic_structure(request):
    """
    Atomic Structure:
    - Bohr energy levels:  Eₙ = -13.6 Z² / n²  (eV)
    - De Broglie wavelength: λ = h / (m v)
    """
    result = None
    error = None
    h = 6.626e-34  # J·s

    if request.method == "POST":
        calc = request.POST.get("calc")
        try:
            if calc == "bohr_energy":
                Z = float(request.POST["Z"])
                n = float(request.POST["n"])
                if n <= 0:
                    error = "Principal quantum number n must be > 0."
                else:
                    E = -13.6 * (Z ** 2) / (n ** 2)
                    result = f"Energy of electron Eₙ = {E:.4f} eV"

            elif calc == "de_broglie":
                m = float(request.POST["m"])
                v = float(request.POST["v"])
                if m <= 0 or v <= 0:
                    error = "Mass and velocity must both be > 0."
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
    Chemical Bonding & Molecular Structure:
    - Approx. % ionic character using Pauling relation:
      %IC ≈ 16(Δχ) + 3.5(Δχ)²,  Δχ = |χA - χB|
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
    - Molarity:      M = n / V
    - Mass percent:  w% = (mass solute / total mass) × 100
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
    Thermodynamics (Chemistry):
    - First law:         ΔU = Q - W
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


# ===================== CHEMISTRY — Syllabus topic pages (placeholders) =====================

def _chem_topic_page(request, title, description):
    """
    Reuse common topic_placeholder for Chemistry topics.
    Design remains same as Physics/Maths placeholder pages.
    """
    return topic_placeholder(
        request,
        subject="chemistry",
        topic=title,
        description=description,
    )


# ---- Physical & General Chemistry ----

def chem_stoichiometry(request):
    return _chem_topic_page(
        request,
        "Stoichiometry & Mole Concept",
        "Mole concept, Avogadro number, empirical and molecular formula, "
        "stoichiometric calculations with mass, volume and moles, limiting reagent and percentage yield."
    )


def chem_states_matter(request):
    return _chem_topic_page(
        request,
        "States of Matter (Gas, Liquid & Solid)",
        "Ideal and real gases, gas laws, kinetic theory basics, intermolecular forces, "
        "liquids and solids – properties and simple numerical ideas."
    )


def chem_periodic_table(request):
    return _chem_topic_page(
        request,
        "Periodic Table & Periodicity",
        "Modern periodic table, classification of elements, trends in atomic/ionic size, "
        "ionisation enthalpy, electron gain enthalpy, electronegativity and valency."
    )


def chem_equilibrium(request):
    return _chem_topic_page(
        request,
        "Chemical & Ionic Equilibrium",
        "Dynamic equilibrium, law of mass action, equilibrium constant, Le Chatelier principle, "
        "pH, buffer solutions, solubility product and common ion effect."
    )


def chem_redox(request):
    return _chem_topic_page(
        request,
        "Redox Reactions",
        "Oxidation–reduction concepts, oxidation number, balancing redox equations, "
        "redox titrations and applications in chemistry."
    )


def chem_electrochemistry(request):
    return _chem_topic_page(
        request,
        "Electrochemistry",
        "Electrolytic and galvanic cells, EMF of cell, Nernst equation, standard electrode potential, "
        "electrochemical series and simple batteries."
    )


def chem_chemical_kinetics(request):
    return _chem_topic_page(
        request,
        "Chemical Kinetics",
        "Rate of reaction, rate law and order, integrated rate expressions for simple reactions, "
        "Arrhenius equation and factors affecting rate."
    )


def chem_surface_chemistry(request):
    return _chem_topic_page(
        request,
        "Surface Chemistry",
        "Adsorption, catalysis, colloids, emulsions and applications of surface phenomena in daily life."
    )


# ---- Inorganic Chemistry ----

def chem_hydrogen(request):
    return _chem_topic_page(
        request,
        "Hydrogen & its Compounds",
        "Position of hydrogen, isotopes, physical and chemical properties, "
        "hydrides, water, hydrogen peroxide and their uses."
    )


def chem_s_block(request):
    return _chem_topic_page(
        request,
        "s-Block Elements",
        "Group 1 and 2 elements – electronic configuration, trends, important compounds "
        "like NaCl, NaOH, CaO, CaCO₃ and their applications."
    )


def chem_p_block(request):
    return _chem_topic_page(
        request,
        "p-Block Elements",
        "Group 13–18 overview, important trends, and main compounds of boron, carbon, nitrogen, "
        "oxygen, halogens and noble gases."
    )


def chem_d_f_block(request):
    return _chem_topic_page(
        request,
        "d- & f-Block Elements",
        "Transition and inner-transition elements, variable oxidation states, colors, "
        "magnetic properties and important alloys."
    )


def chem_coordination(request):
    return _chem_topic_page(
        request,
        "Coordination Compounds",
        "Ligands, coordination number, Werner theory, nomenclature, isomerism and basic applications "
        "like cisplatin, hemoglobin, etc."
    )


def chem_metallurgy(request):
    return _chem_topic_page(
        request,
        "Metallurgy & Qualitative Analysis",
        "General principles of extraction of metals, concentration of ores, "
        "reduction, refining and a short overview of qualitative salt analysis."
    )


def chem_environmental(request):
    return _chem_topic_page(
        request,
        "Environmental Chemistry",
        "Pollution of air, water and soil, ozone depletion, greenhouse effect, "
        "acid rain and strategies for control of pollution."
    )


# ---- Organic Chemistry ----

def chem_organic_basics(request):
    return _chem_topic_page(
        request,
        "Basic Concepts of Organic & IUPAC",
        "Classification of organic compounds, types of bonding and isomerism, "
        "IUPAC nomenclature and basics of reaction mechanisms."
    )


def chem_hydrocarbons(request):
    return _chem_topic_page(
        request,
        "Hydrocarbons (Alkane, Alkene, Alkyne)",
        "Preparation, properties and reactions of alkanes, alkenes and alkynes, "
        "including important named reactions and uses."
    )


def chem_haloalkanes(request):
    return _chem_topic_page(
        request,
        "Haloalkanes & Haloarenes",
        "Nomenclature, nature of C–X bond, substitution and elimination reactions, "
        "and environmental effects of halo-compounds."
    )


def chem_alcohols(request):
    return _chem_topic_page(
        request,
        "Alcohols, Phenols & Ethers",
        "Preparation and properties of alcohols, phenols and ethers with focus on "
        "important reactions and everyday applications."
    )


def chem_aldehydes(request):
    return _chem_topic_page(
        request,
        "Aldehydes, Ketones & Carboxylic Acids",
        "Structure, preparation and characteristic reactions including nucleophilic addition "
        "and acidity of carboxylic acids."
    )


def chem_amines(request):
    return _chem_topic_page(
        request,
        "Amines",
        "Classification, basicity, preparation and important reactions of aliphatic and aromatic amines."
    )


def chem_biomolecules(request):
    return _chem_topic_page(
        request,
        "Biomolecules",
        "Carbohydrates, proteins, vitamins and nucleic acids – basic structures and biological roles."
    )


def chem_polymers(request):
    return _chem_topic_page(
        request,
        "Polymers",
        "Types of polymers (addition, condensation, natural, synthetic), "
        "common examples and their uses in daily life."
    )


def chem_everyday_life(request):
    return _chem_topic_page(
        request,
        "Chemistry in Everyday Life",
        "Drugs, soaps and detergents, food additives, and other chemical substances used in daily life."
    )
