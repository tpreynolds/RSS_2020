function eul = dcm_2_eul(C_I2B)

pitch = asin(-C_I2B(1,3));
roll  = atan2( C_I2B(2,3), C_I2B(3,3) );
yaw   = atan2( C_I2B(1,2), C_I2B(1,1) );

eul = [ roll; pitch; yaw ];

end

