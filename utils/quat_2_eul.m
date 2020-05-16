function eul = quat_2_eul(quat)

dcm = quat_2_dcm_last(quat);
eul = dcm_2_eul(dcm);

end

