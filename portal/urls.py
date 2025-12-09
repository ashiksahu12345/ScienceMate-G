from django.contrib import admin
from django.urls import path, include  
from . import views

urlpatterns = [
    # HOME
    path('', views.home, name='home'),

    # ---------- PHYSICS: Full syllabus topics (placeholders for now) ----------
    path('physics/motion/',             views.physics_motion,          name='physics_motion'),
    path('physics/forces-laws/',        views.physics_forces,          name='physics_forces'),
    path('physics/work-energy-power/',  views.physics_work_energy,     name='physics_work_energy'),
    path('physics/momentum-collisions/',views.physics_momentum,        name='physics_momentum'),
    path('physics/circular-motion/',    views.physics_circular,        name='physics_circular'),
    path('physics/projectile-motion/',  views.physics_projectile,      name='physics_projectile'),
    path('physics/gravitation/',        views.physics_gravity,         name='physics_gravity'),
    path('physics/fluid-mechanics/',    views.physics_fluids,          name='physics_fluids'),
    path('physics/heat-thermodynamics/',views.physics_heat,            name='physics_heat'),
    path('physics/waves-sound/',        views.physics_waves,           name='physics_waves'),
    path('physics/shm-oscillations/',   views.physics_shm,             name='physics_shm'),
    path('physics/ray-optics/',         views.physics_ray_optics,      name='physics_ray_optics'),
    path('physics/wave-optics/',        views.physics_wave_optics,     name='physics_wave_optics'),
    path('physics/electrostatics/',     views.physics_electrostatics,  name='physics_electrostatics'),
    path('physics/electricity-dc/',     views.physics_electricity,     name='physics_electricity'),
    path('physics/magnetism-emi/',      views.physics_magnetism,       name='physics_magnetism'),
    path('physics/ac-circuits/',        views.physics_ac,              name='physics_ac'),
    path('physics/em-waves/',           views.physics_em_waves,        name='physics_em_waves'),
    path('physics/modern-physics/',     views.physics_modern,          name='physics_modern'),
    path('physics/units-dimensions/',   views.physics_units,           name='physics_units'),
    path("physics/units-dimensions-errors/", views.physics_units, name="physics_units"),
    path("physics/ray-optics/", views.physics_ray_optics, name="physics_ray_optics"),
    path("physics/wave-optics/", views.physics_wave_optics, name="physics_wave_optics"),
    path('physics/rotational-motion/',        views.physics_rotational,     name='physics_rotational'),
    path('physics/properties-of-solids/',     views.physics_solids,         name='physics_solids'),
    path('physics/kinetic-theory/',           views.physics_kinetic_theory, name='physics_kinetic_theory'),
    path('physics/semiconductors-electronics/', views.physics_semiconductors, name='physics_semiconductors'),
    path('physics/communication-systems/',    views.physics_communication,  name='physics_communication'),
    path("physics/thermal-properties-of-matter/",      views.physics_thermal_properties,  name="physics_thermal_properties"),
    path("physics/electric-potential-capacitance/",    views.physics_capacitance,         name="physics_capacitance"),
    path("physics/magnetic-effects-of-current/",       views.physics_magnetic_effects,    name="physics_magnetic_effects"),
    path("physics/electromagnetic-induction-detailed/",views.physics_emi_detailed,        name="physics_emi_detailed"),
    path("physics/ac-power-resonance/",                views.physics_ac_power,            name="physics_ac_power"),
    path("physics/laser-physics/",                     views.physics_lasers,              name="physics_lasers"),
    path("physics/optical-fiber-waveguides/",          views.physics_optical_fiber,       name="physics_optical_fiber"),
    path("physics/nuclear-physics-radioactivity/",     views.physics_nuclear,             name="physics_nuclear"),
    path("physics/error-data-analysis/",               views.physics_errors_data,         name="physics_errors_data"),
    path("physics/instruments-measurements/",          views.physics_instruments,         name="physics_instruments"),


# ===================== MATHS ROUTES =====================

# ---- Full syllabus topic pages (main topic screens) ----
path("maths/ap/",                  views.maths_ap,                  name="maths_ap"),
path("maths/real-numbers/",        views.maths_real_numbers,        name="maths_real_numbers"),
path("maths/polynomials/",         views.maths_polynomials,         name="maths_polynomials"),
path("maths/linear-equations/",    views.maths_linear_equations,    name="maths_linear_equations"),
path("maths/quadratic-equations/", views.maths_quadratic_equations, name="maths_quadratic_equations"),
path("maths/ap-topic/",            views.maths_ap,                  name="maths_ap_topic"),  # अगर चाहो तो इसी को इस्तेमाल करो topic button के लिए

# Coordinate & Trigonometry topics
path("maths/coordinate-basics/",   views.maths_coordinate_basics,   name="maths_coordinate_basics"),
path("maths/trig-basic/",          views.maths_trig_basic,          name="maths_trig_basic"),
path("maths/trig-applications/",   views.maths_trig_applications,   name="maths_trig_applications"),

# Geometry of circle & mensuration
path("maths/circles/",             views.maths_circles,             name="maths_circles"),
path("maths/areas-circles/",       views.maths_area_circles,        name="maths_area_circles"),
path("maths/surface-volume/",      views.maths_surface_volume,      name="maths_surface_volume"),

# Algebra of functions & matrices
path("maths/relations-functions/", views.maths_relations_functions, name="maths_relations_functions"),
path("maths/matrices/",            views.maths_matrices,            name="maths_matrices"),
path("maths/determinants/",        views.maths_determinants,        name="maths_determinants"),

# Calculus topics
path("maths/limits-continuity/",   views.maths_limits_continuity,   name="maths_limits_continuity"),
path("maths/differentiation/",     views.maths_differentiation,     name="maths_differentiation"),
path("maths/app-derivatives/",     views.maths_app_derivatives,     name="maths_app_derivatives"),
path("maths/integration/",         views.maths_integration,         name="maths_integration"),
path("maths/app-integrals/",       views.maths_app_integrals,       name="maths_app_integrals"),
path("maths/diff-eq/",             views.maths_diff_eq,             name="maths_diff_eq"),

# 3D, vectors, LPP
path("maths/3d-geometry/",         views.maths_3d_geometry,         name="maths_3d_geometry"),
path("maths/linear-programming/",  views.maths_linear_programming,  name="maths_linear_programming"),
path("maths/vectors/",             views.maths_vectors,             name="maths_vectors"),

# ---- Interactive calculator tools (extra helpers) ----
path("maths/arithmetic/",          views.maths_arithmetic,          name="maths_arithmetic"),
path("maths/algebra/",             views.maths_algebra,             name="maths_algebra"),
path("maths/geometry/",            views.maths_geometry,            name="maths_geometry"),
path("maths/statistics/",          views.maths_statistics,          name="maths_statistics"),
path("maths/probability/",         views.maths_probability,         name="maths_probability"),
path("maths/sphere/",              views.maths_sphere,              name="maths_sphere"),
path("maths/coordinate-midpoint/", views.maths_midpoint,            name="maths_midpoint"),



        # ---------- CHEMISTRY: Interactive + full syllabus topic pages ----------

    # Existing interactive calculators
    path('chemistry/atomic-structure/',      views.chem_atomic_structure,   name='chem_atomic_structure'),
    path('chemistry/chemical-bonding/',      views.chem_chemical_bonding,   name='chem_chemical_bonding'),
    path('chemistry/solutions/',             views.chem_solutions,          name='chem_solutions'),
    path('chemistry/thermodynamics/',        views.chem_thermodynamics,     name='chem_thermodynamics'),

    # Physical & General Chemistry (placeholders for now – same design as physics maths)
    path('chemistry/stoichiometry-mole-concept/',     views.chem_stoichiometry,       name='chem_stoichiometry'),
    path('chemistry/states-of-matter/',               views.chem_states_matter,       name='chem_states_matter'),
    path('chemistry/periodic-table-periodicity/',     views.chem_periodic_table,      name='chem_periodic_table'),
    path('chemistry/chemical-ionic-equilibrium/',     views.chem_equilibrium,         name='chem_equilibrium'),
    path('chemistry/redox-reactions/',                views.chem_redox,               name='chem_redox'),
    path('chemistry/electrochemistry/',               views.chem_electrochemistry,    name='chem_electrochemistry'),
    path('chemistry/chemical-kinetics/',              views.chem_chemical_kinetics,   name='chem_chemical_kinetics'),
    path('chemistry/surface-chemistry/',              views.chem_surface_chemistry,   name='chem_surface_chemistry'),

    # Inorganic Chemistry
    path('chemistry/hydrogen-compounds/',             views.chem_hydrogen,            name='chem_hydrogen'),
    path('chemistry/s-block-elements/',               views.chem_s_block,             name='chem_s_block'),
    path('chemistry/p-block-elements/',               views.chem_p_block,             name='chem_p_block'),
    path('chemistry/d-f-block-elements/',             views.chem_d_f_block,           name='chem_d_f_block'),
    path('chemistry/coordination-compounds/',         views.chem_coordination,        name='chem_coordination'),
    path('chemistry/metallurgy-qualitative-analysis/',views.chem_metallurgy,          name='chem_metallurgy'),
    path('chemistry/environmental-chemistry/',        views.chem_environmental,       name='chem_environmental'),

    # Organic Chemistry
    path('chemistry/organic-basics-iupac/',           views.chem_organic_basics,      name='chem_organic_basics'),
    path('chemistry/hydrocarbons/',                   views.chem_hydrocarbons,        name='chem_hydrocarbons'),
    path('chemistry/haloalkanes-haloarenes/',         views.chem_haloalkanes,         name='chem_haloalkanes'),
    path('chemistry/alcohols-phenols-ethers/',        views.chem_alcohols,            name='chem_alcohols'),
    path('chemistry/aldehydes-ketones-carboxylic-acids/', views.chem_aldehydes,      name='chem_aldehydes'),
    path('chemistry/amines/',                         views.chem_amines,              name='chem_amines'),
    path('chemistry/biomolecules/',                   views.chem_biomolecules,        name='chem_biomolecules'),
    path('chemistry/polymers/',                       views.chem_polymers,            name='chem_polymers'),
    path('chemistry/chemistry-everyday-life/',        views.chem_everyday_life,       name='chem_everyday_life'),



    # ---------- CHATBOT & UPLOAD QUESTION ----------
    path('chatbot/',            views.chatbot,          name='chatbot'),
    path('upload-question/',    views.upload_question,  name='upload_question'),
]
