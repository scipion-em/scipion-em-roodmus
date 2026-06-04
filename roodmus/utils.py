# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     David Herreros (dherreros@cnb.csic.es)
# *
# * National Centre for Biotechnology (CSIC)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


import numpy as np


def normalize_image(image, radius_fraction=0.5):
    image_size = image.shape[0]
    radius_pixels = radius_fraction * image_size

    coords = np.arange(image_size) - (image_size - 1) / 2.0
    xx, yy = np.meshgrid(coords, coords, indexing='ij')
    distance_from_center = np.sqrt(xx ** 2 + yy ** 2)
    outside_circle_mask = distance_from_center > radius_pixels

    # Count the number of valid pixels from the mask (as a float).
    num_noise_pixels = np.sum(outside_circle_mask)

    # Compute mean
    masked_image = np.where(outside_circle_mask, image, 0.0)
    image_sum = np.sum(masked_image)
    mean_val = image_sum / num_noise_pixels

    # Compute standard deviation
    image_sq_sum = np.sum(np.where(outside_circle_mask, np.power(image, 2), 0.0))
    mean_of_squares = image_sq_sum / num_noise_pixels
    variance = mean_of_squares - np.power(mean_val, 2)
    std_val = np.sqrt(np.maximum(0.0, variance))

    # Normalize image
    image = (image - mean_val) / std_val

    return image
