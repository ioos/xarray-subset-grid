import os

import numpy as np
import pytest

from xarray_subset_grid.utils import (
    normalize_bbox_x_coords,
    normalize_polygon_x_coords,
    ray_tracing_numpy,
)

# normalize_polygon_x_coords tests.


def get_test_file_dir():
    """
    returns the test file dir path
    """
    test_file_dir = os.path.join(os.path.dirname(__file__), "example_data")
    return test_file_dir


poly1_180 = np.array(
    [
        [-73, 41],
        [-70, 41],
        [-73, 39],
        [-73, 41],
    ]
)
poly1_360 = np.array(
    [
        [287, 41],
        [290, 41],
        [287, 39],
        [287, 41],
    ]
)

poly2_360 = np.array(
    [
        [234, 41],
        [234, 41],
        [250, 39],
        [290, 41],
    ]
)

poly2_180 = np.array(
    [
        [-126, 41],
        [-126, 41],
        [-110, 39],
        [-70, 41],
    ]
)


@pytest.mark.parametrize(
    "lons, poly, norm_poly",
    [
        ([-85, -84, -83, 10], poly1_180, poly1_180),  # x1
        ([60, 45, 85, 70], poly1_180, poly1_180),  # x2
        ([190, 200, 220, 250, 260], poly1_180, poly1_360),  # x3
        ([-85, -84, -83, 10], poly2_360, poly2_180),  # x1
        ([60, 45, 85, 70], poly2_360, poly2_360),  # x2
        ([190, 200, 220, 250, 260], poly2_360, poly2_360),  # x3
    ],
)
def test_normalize_x_coords(lons, poly, norm_poly):
    lons = np.array(lons)
    normalized_polygon = normalize_polygon_x_coords(lons, np.array(poly))
    print(f"{lons=}")
    print(f"{poly=}")
    print(f"{norm_poly=}")
    print(f"{normalized_polygon=}")

    assert np.allclose(normalized_polygon, norm_poly)


bbox1_180 = [-73, 39, -70, 41]
bbox1_360 = [287, 39, 290, 41]
bbox2_360 = [234, 39, 290, 41]
bbox2_180 = [-126, 39, -70, 41]


@pytest.mark.parametrize(
    "lons, bbox, norm_bbox",
    [
        ([-85, -84, -83, 10], bbox1_180, bbox1_180),  # x1
        ([60, 45, 85, 70], bbox1_180, bbox1_180),  # x2
        ([190, 200, 220, 250, 260], bbox1_180, bbox1_360),  # x3
        ([-85, -84, -83, 10], bbox2_360, bbox2_180),  # x1
        ([60, 45, 85, 70], bbox2_360, bbox2_360),  # x2
        ([190, 200, 220, 250, 260], bbox2_360, bbox2_360),  # x3
    ],
)
def test_normalize_x_coords_bbox(lons, bbox, norm_bbox):
    lons = np.array(lons)
    normalized_polygon = normalize_bbox_x_coords(lons, bbox)
    assert np.allclose(normalized_polygon, norm_bbox)


def test_ray_tracing_numpy():
    """
    minimal test, but at least it'll show it's not totally broken

    NOTE: this function was compared to shapely and a Cython implementation
    """
    poly = [
        (3.0, 3.0),
        (5.0, 8.0),
        (10.0, 5.0),
        (7.0, 1.0),
    ]

    points = np.array(
        [
            (3.0, 6.0),  # outside
            (6.0, 4.0),  # inside
            (9.0, 7.0),  # outside
        ]
    )

    result = ray_tracing_numpy(points[:, 0], points[:, 1], poly)

    assert np.array_equal(result, [False, True, False])
