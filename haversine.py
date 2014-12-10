from math import *

__author__ = 'Steve Kochaver'

def distance(lat_1, lat_2, lon_1, lon_2, r=6367):
    '''
    :param r: Spheriod radius in kilometers. (default = 6367  a pretty good approximation of earth)
    '''

    # Convert decimal degrees to radians
    lat_1, lat_2, lon_1, lon_2 = map(radians, [lat_1, lat_2, lon_1, lon_2])

    # Calculate differences between our latitude and longitude pairs
    dif_lat = abs(lat_1 - lat_2)
    dif_lon = abs(lon_1 - lon_2)

    # Calculate the haversine of our latitiude and longitude differences
    hsin_lat = sin(dif_lat/2)**2
    hsin_lon = sin(dif_lon/2)**2

    # The guts of our lat, lon haversine calculations
    h = hsin_lat + cos(lat_1)*cos(lat_2) * hsin_lon

    # The inverse of the haversine formula
    inv_hsine = 2.0 * asin(sqrt(h))

    # Finally multiply radius by our final inverse haversine calculation
    d = r * inv_hsine

    return d






