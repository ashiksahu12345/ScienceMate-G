from django.urls import path
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

    # ---------- MATHS (already interactive from Part 3A) ----------
        # ---------- MATHS (already interactive from Part 3A) ----------
    path('maths/arithmetic/',           views.maths_arithmetic,        name='maths_arithmetic'),
    path('maths/algebra/',              views.maths_algebra,           name='maths_algebra'),
    path('maths/geometry/',             views.maths_geometry,          name='maths_geometry'),
    path('maths/statistics/',           views.maths_statistics,        name='maths_statistics'),
    path('maths/probability/',          views.maths_probability,       name='maths_probability'),
    path('maths/geometry/sphere/',      views.maths_sphere,            name='maths_sphere'),
    path('maths/coordinate/midpoint/',  views.maths_midpoint,          name='maths_midpoint'),
        # ---------- MATHS: Full syllabus topic pages ----------
    path('maths/real-numbers/',                 views.maths_real_numbers,          name='maths_real_numbers'),
    path('maths/polynomials/',                  views.maths_polynomials,           name='maths_polynomials'),
    path('maths/linear-equations/',             views.maths_linear_equations,      name='maths_linear_equations'),
    path('maths/quadratic-equations/',          views.maths_quadratic_equations,   name='maths_quadratic_equations'),
    path('maths/arithmetic-progressions/',      views.maths_ap,                    name='maths_ap'),

    path('maths/coordinate-basics/',            views.maths_coordinate_basics,     name='maths_coordinate_basics'),
    path('maths/trigonometry-basic/',           views.maths_trig_basic,            name='maths_trig_basic'),
    path('maths/trigonometry-applications/',    views.maths_trig_applications,     name='maths_trig_applications'),
    path('maths/circles/',                      views.maths_circles,               name='maths_circles'),
    path('maths/areas-related-to-circles/',     views.maths_area_circles,          name='maths_area_circles'),
    path('maths/surface-area-volume/',          views.maths_surface_volume,        name='maths_surface_volume'),

    path('maths/relations-functions/',          views.maths_relations_functions,   name='maths_relations_functions'),
    path('maths/matrices/',                     views.maths_matrices,              name='maths_matrices'),
    path('maths/determinants/',                 views.maths_determinants,          name='maths_determinants'),
    path('maths/limits-continuity/',            views.maths_limits_continuity,     name='maths_limits_continuity'),
    path('maths/differentiation/',              views.maths_differentiation,       name='maths_differentiation'),
    path('maths/applications-of-derivatives/',  views.maths_app_derivatives,       name='maths_app_derivatives'),
    path('maths/integration/',                  views.maths_integration,           name='maths_integration'),
    path('maths/applications-of-integrals/',    views.maths_app_integrals,         name='maths_app_integrals'),
    path('maths/differential-equations/',       views.maths_diff_eq,               name='maths_diff_eq'),
    path('maths/vector-algebra/',               views.maths_vectors,               name='maths_vectors'),
    path('maths/3d-geometry/',                  views.maths_3d_geometry,           name='maths_3d_geometry'),
    path('maths/linear-programming/',           views.maths_linear_programming,    name='maths_linear_programming'),
    


    # ---------- CHEMISTRY (already interactive) ----------
    path('chemistry/atomic-structure/',  views.chem_atomic_structure,   name='chem_atomic_structure'),
    path('chemistry/chemical-bonding/',  views.chem_chemical_bonding,   name='chem_chemical_bonding'),
    path('chemistry/solutions/',         views.chem_solutions,          name='chem_solutions'),
    path('chemistry/thermodynamics/',    views.chem_thermodynamics,     name='chem_thermodynamics'),

    # ---------- CHATBOT & UPLOAD QUESTION ----------
    path('chatbot/',            views.chatbot,          name='chatbot'),
    path('upload-question/',    views.upload_question,  name='upload_question'),
]
