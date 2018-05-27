## Shielded applicator from the Lymperopoulou et al. (2004)
##
## Lymperopoulou applicator 
##

:start geometry definition:

    :start geometry:
        library = egs_planes
        name = planes_a
        type = EGS_Zplanes
        positions = -7.0 7.0 9.7
    :stop geometry:

    :start geometry:
        library = egs_cylinders
        name = body_a
        midpoint = 0 0 0
        type = EGS_ZCylinders
        radii = 1.5

        :start media input:
            media = PSU1000
        :stop media input:
    :stop geometry:

    :start geometry:
        library = egs_spheres
        name = end_cap_a
        midpoint = 0 0 7.0
        type = EGS_cSpheres
        radii = 0.8 1.5

        :start media input:
            media = AIR_TG43 PSU1000 
            set medium = 0 0
            set medium = 1 1
        :stop media input:
    :stop geometry:

    :start geometry:
        library = egs_genvelope
        name = full_body_with_hollow_cap
        base geometry = body_a
        incribed geometries = end_cap_a
    :stop geometry:

    :start geometry:
        library = egs_cylinders
        name = rho_coordinates_a
        type = EGS_ZCylinders
        radii = 0.8
    :stop geometry:

    :start geometry:
       library   = egs_planes
       type      = EGS_Zplanes
       positions = -7.1 8.5
       name      = z_coordinates_a
    :stop geometry:

    :start geometry:
       library   = egs_iplanes
       axis      = 0 0 0   0 0 1
       angles    = 90
       name      = phi_coordinates_a
     :stop geometry:

    :start geometry:
        library   = egs_ndgeometry
        name      = shield_a
        dimensions = rho_coordinates_a z_coordinates_a phi_coordinates_a
        :start media input:
           media = AIR_TG43 
        :stop media input:
    :stop geometry:

    :start geometry:
        library = egs_spheres
        name = shield_cap_a
        midpoint = 0 0 7.0
        radii = 0.8
    :stop geometry:

    :start geometry:
        library = egs_iplanes
        axis = 0 0 0   0 0 1 
        angles = 90
        name = hemis_planes_cap_coords
    :stop geometry:

    :start geometry:
        library = egs_ndgeometry
        name = shield_cap
        dimensions = shield_cap_a hemis_planes_cap_coords 

        :start media input:
            media = AIR_TG43 
        :stop media input:
    :stop geometry:

    :start geometry:
        library = egs_cylinders
        name = probe_a
        midpoint = 0 0 0
        type = EGS_ZCylinders
        radii = 0.104 0.16

        :start media input:
            media = AIR_TG43 SS_AISI316L_p8.02
            set medium = 0 0
            set medium = 1 1
        :stop media input:
    :stop geometry:

    :start geometry:
        library = egs_genvelope
        name = shield_and_probe_a
        base geometry = shield_a
        inscribed geometries = probe_a 
    :stop geometry:

    :start geometry:
        library = egs_genvelope
        name = shield_b
        base geometry = shield_cap
        inscribed geometries = probe_a
    :stop geometry:

    :start geometry:
        library = egs_genvelope
        name = shield_c
        base geometry = end_cap_a 
        inscribed geometries = shield_b
    :stop geometry:

    :start geometry:
        library = egs_genvelope
        name = body_w_shield_and_probe_a
        base geometry = full_body_with_hollow_cap
        inscribed geometries = shield_and_probe_a
    :stop geometry:

    :start geometry:
        library = egs_cdgeometry
        name = applicator
        base geometry = planes_a
        set geometry = 0 body_w_shield_and_probe_a
        set geometry = 1 shield_c
    :stop geometry:

    simulation geometry = applicator

:stop geometry definition:

