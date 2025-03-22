# Import global modules. #
import numpy as np
# Import local modules. #
import properties
from parameters_and_constants import num_bins


def interpolate(x1, y1, x2, y2):
    """
    Linearly interpolate the x-coordinate of a given y-coordinate between two points.
    :param x1: x-coordinate of first point.
    :param y1: y-coordinate of first point.
    :param x2: x-coordinate of second point.
    :param y2: y-coordinate of second point.
    :return: Interpolated x-coordinate.
    """

    y = 1.0  # Represents the marginally stable solution with Toomre parameter equal to one.
    return x1 + (y - y1) * (x2 - x1) / (y2 - y1)


def resolve_overlapping_regimes(viscosity, bh_mass, accretion_rate, spin, rotation, x_bins, mask_ges, mask_res,
                                mask_gff):
    """
    If two or more regimes are valid in the same radial bin, assign the bins accordingly to the adjacent regimes.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spin: Dimensionless spin parameter.
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param x_bins: Dimensionless radial bins.
    :param mask_ges: Mask for the 'Gas-ES' radial bins.
    :param mask_res: Mask for the 'Rad-ES' radial bins.
    :param mask_gff: Mask for the 'Gas-FF' radial bins.
    :return:
    """

    # Check if there are overlapping radial bins. If yes, find which regimes overlap and keep the one with the highest
    # pressure. #
    if len(mask_ges) + len(mask_res) + len(mask_gff) > num_bins:
        for overlap_bin in np.intersect1d(mask_ges, mask_res):
            if properties.ratios(x_bins[overlap_bin], viscosity, bh_mass, accretion_rate, spin, rotation, 'Gas-ES')[
                0] < properties.ratios(x_bins[overlap_bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                       'Rad-ES')[0]:

                # Find the index of the overlapping radial bin and remove it. #
                mask_res = np.delete(mask_res, np.where(mask_res == overlap_bin))

            elif properties.ratios(x_bins[overlap_bin], viscosity, bh_mass, accretion_rate, spin, rotation, 'Rad-ES')[
                0] < properties.ratios(x_bins[overlap_bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                       'Gas-ES')[0]:

                # Find the index of the overlapping radial bin and remove it. #
                mask_ges = np.delete(mask_ges, np.where(mask_ges == overlap_bin))

        for overlap_bin in np.intersect1d(mask_ges, mask_gff):
            if properties.ratios(x_bins[overlap_bin], viscosity, bh_mass, accretion_rate, spin, rotation, 'Gas-ES')[
                0] < properties.ratios(x_bins[overlap_bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                       'Gas-FF')[0]:

                # Find the index of the overlapping radial bin and remove it. #
                mask_gff = np.delete(mask_gff, np.where(mask_gff == overlap_bin))

            elif properties.ratios(x_bins[overlap_bin], viscosity, bh_mass, accretion_rate, spin, rotation, 'Gas-FF')[
                0] < properties.ratios(x_bins[overlap_bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                       'Gas-ES')[0]:

                # Find the index of the overlapping radial bin and remove it. #
                mask_ges = np.delete(mask_ges, np.where(mask_ges == overlap_bin))

        for overlap_bin in np.intersect1d(mask_res, mask_gff):
            if properties.ratios(x_bins[overlap_bin], viscosity, bh_mass, accretion_rate, spin, rotation, 'Rad-ES')[
                0] < properties.ratios(x_bins[overlap_bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                       'Gas-FF')[0]:
                # Find the index of the overlapping radial bin and remove it. #
                mask_gff = np.delete(mask_gff, np.where(mask_gff == overlap_bin))

            elif properties.ratios(x_bins[overlap_bin], viscosity, bh_mass, accretion_rate, spin, rotation, 'Gas-FF')[
                0] < properties.ratios(x_bins[overlap_bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                       'Rad-ES')[0]:

                # Find the index of the overlapping radial bin and remove it. #
                mask_res = np.delete(mask_res, np.where(mask_res == overlap_bin))

    return mask_ges, mask_res, mask_gff if len(mask_ges) + len(mask_res) + len(mask_gff) <= num_bins \
        else (print('Terminating "resolve_overlapping_regimes"! Overlapping radial bins still exist', sep=''), exit())


def resolve_empty_regimes(viscosity, bh_mass, accretion_rate, spin, rotation, x_bins, mask_ges, mask_res, mask_gff):
    """
    If a radial bin is empty, assign it to the adjacent regime with the lowest 'pressure_ratio'.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spin: Dimensionless spin parameter.
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param x_bins: Dimensionless radial bins.
    :param mask_ges: Mask for the 'Gas-ES' radial bins.
    :param mask_res: Mask for the 'Rad-ES' radial bins.
    :param mask_gff: Mask for the 'Gas-FF' radial bins.
    :return:
    """
    if len(mask_ges) + len(mask_res) + len(mask_gff) < num_bins:
        # Find the index of the empty radial bin(s) by comparing all mask with an array of length 'num_bins'. #
        indices = np.arange(0, num_bins)
        all_masks = np.concatenate((mask_ges, mask_res, mask_gff))
        indices = np.setdiff1d(indices, all_masks)

        for index in indices:
            # Check if the empty bin is at the ISCO. If yes, then assign the ISCO bin to the regime that is valid in
            # the second radial bin. #
            if index == 0:
                if 1 in mask_ges: mask_ges = np.insert(mask_ges, 0, 0)
                if 1 in mask_res: mask_res = np.insert(mask_res, 0, 0)
                if 1 in mask_gff: mask_gff = np.insert(mask_gff, 0, 0)

            # Check if the empty bin is at the outermost edge. If yes, then assign the ISCO bin to the regime that is
            # valid in the second to last radial bin. #
            if index == num_bins - 1:
                if num_bins - 1 in mask_ges: mask_ges = np.insert(mask_ges, num_bins, num_bins - 1)
                if num_bins - 1 in mask_res: mask_res = np.insert(mask_res, num_bins, num_bins - 1)
                if num_bins - 1 in mask_gff: mask_gff = np.insert(mask_gff, num_bins, num_bins - 1)

            # Check if the empty bin is inbetween the ISCO and the outermost edge. If yes, then find the regime that has
            # the minimum 'pressure_ratio' between all valid regimes in the bin just before and right after and assign
            # it to the previously empty bin. #
            if 0 < index < num_bins - 1:
                # Check if the empty bin is inbetween two regions where the same regime is valid. #
                if index - 1 in mask_ges and index + 1 in mask_ges: mask_ges = np.insert(mask_ges, index, index)
                if index - 1 in mask_res and index + 1 in mask_res: mask_res = np.insert(mask_res, index, index)
                if index - 1 in mask_gff and index + 1 in mask_gff: mask_gff = np.insert(mask_gff, index, index)

                # Check if the empty bin is inbetween two different regimes. #
                if (index - 1 in mask_ges and index + 1 in mask_res) \
                        or (index - 1 in mask_res and index + 1 in mask_ges):
                    if properties.ratios(x_bins[index], viscosity, bh_mass, accretion_rate, spin, rotation, 'Gas-ES')[
                        0] < properties.ratios(x_bins[index], viscosity, bh_mass, accretion_rate, spin, rotation,
                                               'Rad-ES')[0]:
                        # Avoid IndexError in case where an index moves back and forth between two regimes. #
                        try:
                            mask_ges = np.insert(mask_ges, index, index)
                        except IndexError:
                            mask_res = np.insert(mask_res, index, index)

                    elif properties.ratios(x_bins[index], viscosity, bh_mass, accretion_rate, spin, rotation, 'Rad-ES')[
                        0] < properties.ratios(x_bins[index], viscosity, bh_mass, accretion_rate, spin, rotation,
                                               'Gas-ES')[0]:
                        # Avoid IndexError in case where an index moves back and forth between two regimes. #
                        try:
                            mask_res = np.insert(mask_res, index, index)
                        except IndexError:
                            mask_ges = np.insert(mask_ges, index, index)

                if (index - 1 in mask_ges and index + 1 in mask_gff) \
                        or (index - 1 in mask_gff and index + 1 in mask_ges):
                    if properties.ratios(x_bins[index], viscosity, bh_mass, accretion_rate, spin, rotation, 'Gas-ES')[
                        0] < properties.ratios(x_bins[index], viscosity, bh_mass, accretion_rate, spin, rotation,
                                               'Gas-FF')[0]:
                        # Avoid IndexError in cases where an index moves back and forth between two regimes. #
                        try:
                            mask_ges = np.insert(mask_ges, index, index)
                        except IndexError:
                            mask_gff = np.insert(mask_gff, index, index)

                    elif properties.ratios(x_bins[index], viscosity, bh_mass, accretion_rate, spin, rotation, 'Gas-FF')[
                        0] < properties.ratios(x_bins[index], viscosity, bh_mass, accretion_rate, spin, rotation,
                                               'Gas-ES')[0]:
                        # Avoid IndexError in case where an index moves back and forth between two regimes. #
                        try:
                            mask_gff = np.insert(mask_gff, index, index)
                        except IndexError:
                            mask_ges = np.insert(mask_ges, index, index)

                if (index - 1 in mask_res and index + 1 in mask_gff) \
                        or (index - 1 in mask_gff and index + 1 in mask_res):
                    if properties.ratios(x_bins[index], viscosity, bh_mass, accretion_rate, spin, rotation, 'Rad-ES')[
                        0] < properties.ratios(x_bins[index], viscosity, bh_mass, accretion_rate, spin, rotation,
                                               'Gas-FF')[0]:
                        # Avoid IndexError in case where an index moves back and forth between two regimes. #
                        try:
                            mask_res = np.insert(mask_res, index, index)
                        except IndexError:
                            mask_gff = np.insert(mask_gff, index, index)

                    elif properties.ratios(x_bins[index], viscosity, bh_mass, accretion_rate, spin, rotation, 'Gas-FF')[
                        0] < properties.ratios(x_bins[index], viscosity, bh_mass, accretion_rate, spin, rotation,
                                               'Rad-ES')[0]:
                        # Avoid IndexError in case where an index moves back and forth between two regimes. #
                        try:
                            mask_gff = np.insert(mask_gff, index, index)
                        except IndexError:
                            mask_res = np.insert(mask_res, index, index)

    return mask_ges, mask_res, mask_gff if len(mask_ges) + len(mask_res) + len(mask_gff) >= num_bins \
        else (print('Terminating "resolve_empty_regimes"! Empty radial bin(s) still exist', sep=''), exit())


def resolve_lone_regimes(viscosity, bh_mass, accretion_rate, spin, rotation, x_bins, mask_ges, mask_res,
                         mask_gff):
    """
    If a regime is valid in only one radial bin, assign the bin accordingly to an adjacent regimes.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spin: Dimensionless spin parameter.
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param x_bins: Dimensionless radial bins.
    :param mask_ges: Mask for the 'Gas-ES' radial bins.
    :param mask_res: Mask for the 'Rad-ES' radial bins.
    :param mask_gff: Mask for the 'Gas-FF' radial bins.
    :return:
    """

    # Check if 'Gas-ES' is anywhere valid.
    if len(mask_ges) > 1:
        # Check if 'Gas-ES' is valid at the ISCO but not in the second bin. If yes, then extend the adjacent valid
        # regime to 'r_isco' to avoid shrinking the accretion disc.
        if mask_ges[0] == 0 and mask_ges[1] != 1:
            if 1 in mask_res:
                mask_ges = np.delete(mask_ges, 0)
                mask_res = np.insert(mask_res, 0, 0)
            elif 1 in mask_gff:
                mask_ges = np.delete(mask_ges, 0)
                mask_gff = np.insert(mask_gff, 0, 0)
        # Check if 'Gas-ES' is valid in the last but not the second to last bin. If yes, then extend the adjacent valid
        # regime to 'r_outermost' to avoid shrinking the accretion disc. #
        if mask_ges[-1] == num_bins - 1 and mask_ges[-2] != num_bins - 2:
            if num_bins - 1 in mask_res:
                mask_ges = np.delete(mask_ges, -1)
                mask_res = np.insert(mask_res, num_bins, num_bins - 1)
            elif num_bins - 1 in mask_gff:
                mask_ges = np.delete(mask_ges, -1)
                mask_gff = np.insert(mask_gff, num_bins, num_bins - 1)

    # Check if 'Rad-ES' is anywhere valid.
    if len(mask_res) > 1:
        # Check if 'Rad-ES' is valid at the ISCO but not in the second bin. If yes, then extend the adjacent valid
        # regime to 'r_isco' to avoid shrinking the accretion disc.
        if mask_res[0] == 0 and mask_res[1] != 1:
            if 1 in mask_ges:
                mask_res = np.delete(mask_res, 0)
                mask_ges = np.insert(mask_ges, 0, 0)
            elif 1 in mask_gff:
                mask_res = np.delete(mask_res, 0)
                mask_gff = np.insert(mask_gff, 0, 0)
        # Check if 'Rad-ES' is valid in the last but not the second to last bin. If yes, then extend the adjacent valid
        # regime to 'r_outermost' to avoid shrinking the accretion disc. #
        if mask_res[-1] == num_bins - 1 and mask_res[-2] != num_bins - 2:
            if num_bins - 1 in mask_ges:
                mask_res = np.delete(mask_res, -1)
                mask_ges = np.insert(mask_ges, num_bins, num_bins - 1)
            elif num_bins - 1 in mask_gff:
                mask_res = np.delete(mask_res, -1)
                mask_gff = np.insert(mask_gff, num_bins, num_bins - 1)

    # Check if 'Gas-FF' is anywhere valid.
    if len(mask_gff) > 1:
        # Check if 'Gas-FF' is valid at the ISCO but not in the second bin. If yes, then extend the adjacent valid
        # regime to 'r_isco' to avoid shrinking the accretion disc.
        if mask_gff[0] == 0 and mask_gff[1] != 1:
            if 1 in mask_ges:
                mask_gff = np.delete(mask_gff, 0)
                mask_ges = np.insert(mask_ges, 0, 0)
            elif 1 in mask_res:
                mask_gff = np.delete(mask_gff, 0)
                mask_res = np.insert(mask_res, 0, 0)
        # Check if 'Gas-FF' is valid in the last but not the second to last bin. If yes, then extend the adjacent valid
        # regime to 'r_outermost' to avoid shrinking the accretion disc. #
        if mask_gff[-1] == num_bins - 1 and mask_gff[-2] != num_bins - 2:
            if num_bins - 1 in mask_ges:
                mask_gff = np.delete(mask_gff, -1)
                mask_ges = np.insert(mask_ges, num_bins, num_bins - 1)
            elif num_bins - 1 in mask_res:
                mask_gff = np.delete(mask_gff, -1)
                mask_res = np.insert(mask_res, num_bins, num_bins - 1)

    return mask_ges, mask_res, mask_gff
