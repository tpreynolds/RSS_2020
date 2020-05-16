function [C_I2B,dC_deul_I2B,dC_deul_B2I] = eul_2_dcm(eul)

roll = eul(1);
pitch = eul(2);
yaw = eul(3);

cr = cos(roll);
sr = sin(roll);
cp = cos(pitch);
sp = sin(pitch);
cy = cos(yaw);
sy = sin(yaw);

C_I2B = [ cp*cy,            cp*sy,              -sp;
          sr*sp*cy-cr*sy,   sr*sp*sy+cr*cy,     sr*cp;
          cr*sp*cy+sr*sy,   cr*sp*sy-sr*cy,     cr*cp ];
      
dC_deul_I2B = zeros(3,3,3);

dC_deul_I2B(:,:,1) = [ 0.0,                 -sp*cy,     -cp*sy;
                       cr*sp*cy+sr*sy,      sr*cp*cy,   -sr*sp*sy-cr*cy;
                       -sr*sp*cy+cr*sy,     cr*cp*cy,   -cr*sp*sy+sr*cy ];
                   
dC_deul_I2B(:,:,2) = [ 0.0,                 -sp*sy,     cp*cy;
                       cr*sp*sy-sr*cy,      sr*cp*sy,   sr*sp*cy-cr*sy;
                       -sr*sp*sy-cr*cy,     cr*cp*sy,   cr*sp*cy+sr*sy ];
                   
dC_deul_I2B(:,:,3) = [ 0.0,     -cp,        0.0;
                       cr*cp,   -sr*sp,     0.0;
                      -sr*cp,   -cr*sp,     0.0 ];

dC_deul_B2I = zeros(3,3,3);

dC_deul_B2I(:,:,1) = [ 0.0,     -sp*cy,     -cp*sy;
                       0.0,     -sp*sy,     cp*cy;
                       0.0,     -cp,        0.0 ];

dC_deul_B2I(:,:,2) = [ cr*sp*cy+sr*sy,  sr*cp*cy,   -sr*sp*sy-cr*cy;
                       cr*sp*sy-sr*cy,  sr*cp*sy,   sr*sp*cy-cr*sy;
                       cr*cp,           -sr*sp,     0.0 ];
                   
dC_deul_B2I(:,:,3) = [ -sr*sp*cy+cr*sy, cr*cp*cy,   -cr*sp*sy+sr*cy;
                       -sr*sp*sy-cr*cy, cr*cp*sy,   cr*sp*cy+sr*sy;
                       -sr*cp,          -cr*sp,     0.0 ];

end