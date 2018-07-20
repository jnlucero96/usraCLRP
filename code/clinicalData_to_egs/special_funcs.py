#!/usr/bin/env
# author: Joseph Lucero
# created on: 28 June 2018 17:49:24
# purpose: functions that do special tasks


def perform_metallic_artifact_reduction(
    seed_locations, slice_lst_array, BOUNDS, VOXEL_DIMS, VOXEL_CENTERS,
    INTERCEPT, SIZE_OF_GRID
):

    print
    print '-' * 50
    print '10' * 10 + "Attempting to perform Metallic Artifact Reduction"
    print '-' * 50
    print

    # These values are hardcoded but not sure why. Should these be made
    # accessible to the user?
    replace = -90.2
    cutoff = 200
    cyl_max_radius = 0.5  # why is there a cylinder defined?
    mar_counter = 0

    x_bounds, y_bounds, z_bounds = BOUNDS

    for seed in seed_locations:
        z_pos = int((seed[2] - (10 * z_bounds[0])) // VOXEL_DIMS[2])
        for index, center in enumerate(VOXEL_CENTERS):
            x_displacement = abs(center[0] - (seed[0] / 10)) - 0.05 * VOXEL_DIMS
            y_displacement = abs(center[0] - (seed[0] / 10)) - 0.05 * VOXEL_DIMS

            if (x_displacement**2 + y_displacement**2) < cyl_max_radius**2:
                x_pos = index % SIZE_OF_GRID[0]
                y_pos = index // SIZE_OF_GRID[0]

                # This indexing is inefficient...
                # Could we possibly turn this into a numpy array without
                # too much trouble?
                if slice_lst_array[z_pos+2][y_pos][x_pos] + INTERCEPT > cutoff:
                    mar_counter += 1
                    slice_lst_array[z_pos+2][y_pos][x_pos] = replace - INTERCEPT

                if slice_lst_array[z_pos+1][y_pos][x_pos] + INTERCEPT > cutoff:
                    mar_counter += 1
                    slice_lst_array[z_pos+1][y_pos][x_pos] = replace - INTERCEPT

                if slice_lst_array[z_pos][y_pos][x_pos] + INTERCEPT > cutoff:
                    mar_counter += 1
                    slice_lst_array[z_pos][y_pos][x_pos] = replace - INTERCEPT

                if slice_lst_array[z_pos+1][y_pos][x_pos] + INTERCEPT > cutoff:
                    mar_counter += 1
                    slice_lst_array[z_pos+1][y_pos][x_pos] = replace - INTERCEPT

                if slice_lst_array[z_pos+2][y_pos][x_pos] + INTERCEPT > cutoff:
                    mar_counter += 1
                    slice_lst_array[z_pos+2][y_pos][x_pos] = replace - INTERCEPT

    print "Metallic Artifact Reduction complete. {0} high density voxels \
    replaced".format(mar_counter)

    return slice_lst_array
