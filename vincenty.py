from math import *

__author__ = 'Steve Kochaver'

def inverse_solution(lat_1, lat_2, lon_1, lon_2, a=6378137.0, f=1.0/298.257223563):
    # Semi-major axis length of our oblate spheroid in meters. 6378137.0 for WGS-84.
    # a = 6378137.0

    # Flattening of ellipsoid. 1/298.257223563 for WGS-84.
    # f = 1.0/298.257223563

    # Semi-minor axis length of the ellipsoid
    b = (1 - f) * a

    # Convert decimal degrees to radians
    lat_1, lat_2, lon_1, lon_2 = map(radians, [lat_1, lat_2, lon_1, lon_2])

    # Calculate reduced latitiudes.
    rlat_1 = atan((1-f)*tan(lat_1))

    rlat_2 = atan((1-f)*tan(lat_2))

    # Converging lambda variable to express our distance error. Starts as
    # Longitude difference.
    L = lon_2 - lon_1

    lambda_conv = L  # First approximation
    lambda_diff = 1  # This will be our test of convergence (arbitrary 1)
    iter_lim = 0  # Iteration limit before declared non-converging

    while lambda_diff > (10**-12):

        # Calculate the sine of the arc length for use in a later calculation
        sin_arc_len = sqrt((cos(rlat_2)*sin(lambda_conv))**2 + (cos(rlat_1)*sin(rlat_2)-sin(rlat_1) *
                                                                cos(rlat_2)*cos(lambda_conv))**2)

        # Break loop if points are co-incident (thanks movable type scripts)
        if sin_arc_len == 0:
            return 0, 0, 0

        # Calculate the cosine of the arc length for use in a later calculation
        cos_arc_len = sin(rlat_1)*sin(rlat_2) + cos(rlat_1)*cos(rlat_2)*cos(lambda_conv)

        # Using the arctangent of the sin and cosine of the arclength calculate the arclength
        arc_len = atan2(sin_arc_len, cos_arc_len)

        # Calculate sine of azimuth of the geodesic at equator (for further trig)
        sin_az_eq = (cos(rlat_1)*cos(rlat_2)*sin(lambda_conv)) / sin_arc_len

        # Calculate the squared cosine of the geodesic equatorial aziumuth. I promise this is leading somewhere.
        cos_sq_az_eq = 1 - sin_az_eq**2

        # Calculate cosine of twice the arc distance on the sphere form the equator to midpoint of connecting line
        try:
            # If we're dividing by zero we're on the equatorial line (thus the try)
            cos_arc_len_eq = cos_arc_len - (2*sin(rlat_1)*sin(rlat_2) / cos_sq_az_eq)

        except:
            # cos_sq_az_eq will equal zero if on equatorial line
            cos_arc_len_eq = 0

        # Calculate a constant Vincenty uses for convergence
        C = (f/16.0) * cos_sq_az_eq * (4 + f * (4 - 3 * cos_sq_az_eq))

        # Finally apply all our hard work to this iteration of the converging lambda variable
        old_lambda_conv = lambda_conv
        lambda_conv = L + (1 - C) * f * sin_az_eq * (arc_len + C * sin_arc_len *
                                                     (cos_arc_len_eq + C*cos_arc_len *
                                                      (-1 + 2 * cos_arc_len_eq**2)))

        lambda_diff = abs(lambda_conv - old_lambda_conv)

        # Here is where we manage our iteration limiting. For most cases no reasonable attempt should go over 250 iterations.
        iter_lim += 1
        if iter_lim > 250:
            raise ValueError('Failed to converge')

    # Caclulate the U squared value constant
    u_sq = cos_sq_az_eq * ((a**2 - b**2) / b**2)

    # Calculate an A constant specific to the this Vincenty calculation
    A = 1.0 + (u_sq / 16384.0) * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0*u_sq)))

    # Calculate B constant. I know there are a lot of these but Vincenty was a smart guy
    B = (u_sq / 1024.0) * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)))

    # Calculate the change in arc length from our extensive calculations
    delta_arc_len = B * sin_arc_len * (cos_arc_len_eq + (B/4.0) *
                                       (cos_arc_len * (-1.0 + 2.0 * cos_arc_len_eq**2) - (B/6.0) * cos_arc_len_eq *
                                        (-3+4 * sin_arc_len**2) * (-3+4*cos_arc_len_eq**2)))

    # Calculate the arc length between the two points (in units of b)
    s = b*A*(arc_len - delta_arc_len)

    # Calculate the forward azimuth of point 1
    fwd_az_1 = atan2(cos(rlat_2)*sin(lambda_conv), (cos(rlat_1)*sin(rlat_2) - sin(rlat_1)*cos(rlat_2)*cos(lambda_conv)))

    # Calculate the forward azimuth of point 2
    fwd_az_2 = atan2(cos(rlat_1)*sin(lambda_conv), (-sin(rlat_1)*cos(rlat_2) + cos(rlat_1)*sin(rlat_2)*cos(lambda_conv)))

    fwd_az_1 = (fwd_az_1 + (2*pi)) % (2*pi)
    fwd_az_2 = (fwd_az_2 + (2*pi)) % (2*pi)

    fwd_az_1 = degrees(fwd_az_1)
    fwd_az_2 = degrees(fwd_az_2)

    return s, fwd_az_1, fwd_az_2

def direct_solution(lat_1, lon_1, s, fwd_az_1, a=6378137.0, f=1.0/298.257223563):

    # Convert latitude longitude the initial forward azimuth into radians
    lat_1, lon_1, fwd_az_1 = map(radians, [lat_1, lon_1, fwd_az_1])

    # Semi-minor axis length of the ellipsoid
    b = (1 - f) * a

    # Calculate reduced latitude for our inital point's latitude
    rlat_1 = atan((1-f)*tan(lat_1))

    # This is the initial arc length variable of our point
    arc_len_1 = atan2(tan(rlat_1), cos(fwd_az_1))

    # Making a variable for the sine of the forward azimuth makes future calculation a little cleaner
    sin_fwd_az = cos(rlat_1)*sin(fwd_az_1)

    # Another variable for the squared cosine of the forward azimuth of the equation
    cos_sq_fwd_az = (1 - sin(fwd_az_1))*(1 + sin(fwd_az_1))

    # Calculate the U squared variable as in the inverse solution
    u_sq = cos_sq_fwd_az * ((a**2 - b**2) / b**2)

    # Calculate an A constant specific to the this Vincenty calculation
    A = 1.0 + (u_sq / 16384.0) * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0*u_sq)))

    # Calculate B constant. Again, specific to the equation.
    B = (u_sq / 1024.0) * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)))

    # Set up the initial variables for our convergence loop
    arc_len = s / (b * A)
    arc_len_diff = 1
    iter_lim = 0

    while arc_len_diff > 10**-12:

        # Calculate twice the arc length m variable.
        two_arc_m = 2*arc_len_1 + arc_len

        # Here we begin calculating the difference for the convergence to be applied.
        delta_arc_len = B * sin(arc_len) * (cos(two_arc_m) + (B/4) * (cos(arc_len) * (-1 + 2 * cos(two_arc_m)**2) -
                                                                      (B/6)*cos(two_arc_m)*(-3+4*sin(arc_len)**2) *
                                                                      (-3+4*cos(two_arc_m)**2)))

        # Store the old arc length, calculate the new, and check the difference until we reach appropriate error in our convergence
        old_arc_len = arc_len
        arc_len = (s / (b * A)) + delta_arc_len

        arc_len_diff = abs(arc_len - old_arc_len)

        # Here is where we manage our iteration limiting. For most cases no reasonable attempt should go over 250 iterations.

        iter_lim += 1
        if iter_lim > 250:
            raise ValueError('Failed to converge')

    # Calculate the latitude of our destination point
    lat_2 = atan2(sin(rlat_1)*cos(arc_len) + cos(rlat_1)*sin(arc_len)*cos(fwd_az_1),
                  (1-f)*sqrt(sin_fwd_az**2 + (sin(rlat_1)*sin(arc_len) - cos(rlat_1) *
                                              cos(arc_len)*cos(fwd_az_1))**2))
    # Calculate the difference variable we will use to find the destination latitude
    lambda_x = atan2(sin(arc_len)*sin(fwd_az_1),
                     cos(rlat_1)*cos(arc_len) - sin(rlat_1) *
                     sin(arc_len)*cos(fwd_az_1))

    # Nice variables and identities that Vincenty had all figured out for this stuff
    C = (f/16)*cos_sq_fwd_az*(4 + f*(4 - 3 * cos_sq_fwd_az))

    L = lambda_x - (1 - C) * f * sin_fwd_az * (arc_len + C * sin(arc_len) *
                                              (cos(two_arc_m) + C * cos(arc_len) *
                                               (-1 + 2 * cos(two_arc_m)**2)))

    # Calculate the destination longitude and put it on a -180 to 180 format
    lon_2 = (lon_1 + L + 3*pi) % (2*pi) - pi

    # Calculate the forward azimuth after reaching the destination point
    fwd_az_2 = atan2(sin_fwd_az, -(sin(rlat_1)*sin(arc_len) - cos(rlat_1)*cos(arc_len)*cos(fwd_az_1)))

    # Cast all our variables back into degrees and put onto a 0 to 360 scale in the case of our freshly calculated azimuth
    fwd_az_2 = degrees(fwd_az_2 + 2*pi) % (2*pi)
    lat_2 = degrees(lat_2)
    lon_2 = degrees(lon_2)

    return lat_2, lon_2, fwd_az_2