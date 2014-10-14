# The code was originally written in JavaScrpt by Chris Veness:
# http://www.movable-type.co.uk/scripts/latlong-convert-coords.html
# http://www.movable-type.co.uk/scripts/latlong-gridref.html
# And was originally realeased on the above site under LGPL license:
# (c) Chris Veness 2005-2010
# But is now released on the above site under CC-BY
# (c) Chris Veness 2005-2012
#
module Breasal
  class WGS84
    include Breasal::Utils

    def initialize(options)
      @latitude  = options[:latitude]
      @longitude = options[:longitude]

      #ellipse parameters
      @e = { :wgs84 =>   { :a=> 6378137,     :b=> 6356752.3142, :f=> 1 / 298.257223563 },
            :airy1830 => { :a=> 6377563.396, :b=> 6356256.910,  :f=> 1 / 299.3249646   } };

      #helmert transform parameters
      @h = { :wgs84toOSGB36 => { :tx=> -446.448,  :ty=>  125.157,  :tz=> -542.060,  # m
                                 :rx=>  -0.1502,  :ry=>   -0.2470, :rz=>  -0.8421,  # sec
                                 :s=>   20.4894 },                                  # ppm
             :osgb36toWGS84 => { :tx=>  446.448,  :ty=> -125.157,  :tz=>  542.060,
                                 :rx=>    0.1502, :ry=>    0.2470, :rz=>    0.8421,
                                 :s=>   -20.4894 } };
    end

    def to_OSGB36
      convert(@latitude, @longitude, 0, @e[:wgs84], @h[:wgs84toOSGB36], @e[:airy1830]);
    end

    # OSGB36 lat lon to OS UK grid eastings & northings
    # http://www.movable-type.co.uk/scripts/latlong-gridref.html
    def to_UKGrid
      (oslat, oslon, _) = self.to_OSGB36

      lat = deg_to_rad(oslat);
      lon = deg_to_rad(oslon);

      a = 6377563.396; b = 6356256.910           # Airy 1830 major & minor semi-axes
      f0 = 0.9996012717                          # NatGrid scale factor on central meridian
      lat0 = deg_to_rad(49);lon0 = deg_to_rad(-2)# NatGrid true origin
      n0 = -100000; e0 = 400000;                 # northing & easting of true origin, metres
      e2 = 1 - (b*b) / (a*a);                    # eccentricity squared
      n = (a-b) / (a+b); n2 = n*n; n3 = n*n*n;

      cosLat = Math.cos(lat); sinLat = Math.sin(lat);
      nu = a*f0/Math.sqrt(1-e2*sinLat*sinLat);              # transverse radius of curvature
      rho = a*f0*(1-e2) / ( (1-e2*sinLat*sinLat) ** 1.5);  # meridional radius of curvature
      eta2 = nu/rho-1;

      ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (lat-lat0);
      mb = (3*n + 3*n*n + (21/8)*n3) * Math.sin(lat-lat0) * Math.cos(lat+lat0);
      mc = ((15/8)*n2 + (15/8)*n3) * Math.sin(2*(lat-lat0)) * Math.cos(2*(lat+lat0));
      md = (35/24)*n3 * Math.sin(3*(lat-lat0)) * Math.cos(3*(lat+lat0));
      m = b * f0 * (ma - mb + mc - md);              # meridional arc

      cos3lat = cosLat*cosLat*cosLat;
      cos5lat = cos3lat*cosLat*cosLat;
      tan2lat = Math.tan(lat)*Math.tan(lat);
      tan4lat = tan2lat*tan2lat;

      i = m + n0;
      ii = (nu/2)*sinLat*cosLat;
      iii = (nu/24)*sinLat*cos3lat*(5-tan2lat+9*eta2);
      iiiA = (nu/720)*sinLat*cos5lat*(61-58*tan2lat+tan4lat);
      iv = nu*cosLat;
      v = (nu/6)*cos3lat*(nu/rho-tan2lat);
      vi = (nu/120) * cos5lat * (5 - 18*tan2lat + tan4lat + 14*eta2 - 58*tan2lat*eta2);

      dLon = lon-lon0;
      dLon2 = dLon*dLon
      dLon3 = dLon2*dLon
      dLon4 = dLon3*dLon
      dLon5 = dLon4*dLon
      dLon6 = dLon5*dLon

      n = i + ii*dLon2 + iii*dLon4 + iiiA*dLon6;
      e = e0 + iv*dLon + v*dLon3 + vi*dLon5;

      return { easting: e.round, northing: n.round }   #return raw easting and northings instead
    end

    private

    def convert(p1lat, p1lon, p1height, e1, t, e2)
       # -- convert polar to cartesian coordinates (using ellipse 1)

       p1lat = deg_to_rad(p1lat); p1lon = deg_to_rad(p1lon);

       a = e1[:a]; b = e1[:b];

       sinPhi = Math.sin(p1lat); cosPhi = Math.cos(p1lat);
       sinLambda = Math.sin(p1lon); cosLambda = Math.cos(p1lon);
       h = p1height;

       eSq = (a*a - b*b) / (a*a);
       nu = a / Math.sqrt(1 - eSq*sinPhi*sinPhi);

       x1 = (nu+h) * cosPhi * cosLambda;
       y1 = (nu+h) * cosPhi * sinLambda;
       z1 = ((1-eSq)*nu + h) * sinPhi;


       # -- apply helmert transform using appropriate params

       tx = t[:tx]; ty = t[:ty]; tz = t[:tz];
       rx = t[:rx] / 3600 * Math::PI/180;  #normalise seconds to radians
       ry = t[:ry] / 3600 * Math::PI/180;
       rz = t[:rz] / 3600 * Math::PI/180;
       s1 = t[:s] / 1e6 + 1;              #normalise ppm to (s+1)

       #apply transform
       x2 = tx + x1*s1 - y1*rz + z1*ry;
       y2 = ty + x1*rz + y1*s1 - z1*rx;
       z2 = tz - x1*ry + y1*rx + z1*s1;


       # -- convert cartesian to polar coordinates (using ellipse 2)

       a = e2[:a]; b = e2[:b];
       precision = 4 / a;  # results accurate to around 4 metres

       eSq = (a*a - b*b) / (a*a);
       p = Math.sqrt(x2*x2 + y2*y2);
       phi = Math.atan2(z2, p*(1-eSq)); phiP = 2 * Math::PI;
       while ( (phi-phiP).abs > precision) do
          nu = a / Math.sqrt(1 - eSq*Math.sin(phi)*Math.sin(phi));
          phiP = phi;
          phi = Math.atan2(z2 + eSq*nu*Math.sin(phi), p);
       end
       lambda = Math.atan2(y2, x2);
       h = p/Math.cos(phi) - nu;

       #return array [lat,lon,height]
       return [ rad_to_deg(phi), rad_to_deg(lambda), h ];
    end
  end
end
