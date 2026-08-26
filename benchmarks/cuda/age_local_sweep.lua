-- Final-system exports for every AGE topology in a narrow continuous bucket.
-- Run in a directory containing base_motor.fem with
-- XFEMM_LINEAR_SYSTEM_EXPORT_BY_TOPOLOGY set.

showconsole()
open("base_motor.fem")

DEG = 3.141592653589793 / 180
Npoles = 12

function set_balanced_delta_iq(Iq_line_peak, theta_mech)
    theta_line = (Npoles/2) * theta_mech - 30
    line_a = Iq_line_peak * sin(theta_line * DEG)
    line_b = Iq_line_peak * sin((theta_line - 120) * DEG)
    line_c = Iq_line_peak * sin((theta_line + 120) * DEG)
    mi_modifycircprop("A", 1, (line_a - line_b) / 3)
    mi_modifycircprop("B", 1, (line_b - line_c) / 3)
    mi_modifycircprop("C", 1, (line_c - line_a) / 3)
end

angles = {0, 0.4, 0.8, 1.2}
for index = 1, 4 do
    angle = angles[index]
    set_balanced_delta_iq(1.0, angle)
    mi_modifyboundprop("SlidingBand", 11, angle)
    print("XFEMM_AGE_LOCAL_BEGIN," .. angle)
    mi_analyze(1)
    print("XFEMM_AGE_LOCAL_END," .. angle)
end

quit()
