-- Postprocess GPU-injected .ans files with XFEMM's normal physics code.
-- Run from a directory containing gpu_age_{0,60,120}.{fem,ans}.

showconsole()

angles = {0, 60, 120}
print("XFEMM_GPU_AGE_HEADER,angle,torque,lambda_a,lambda_b,lambda_c," ..
      "gap_a_0,h1_brc,h1_brs,h1_btc,h1_bts,h6_brc,h6_brs,h6_btc,h6_bts")

for index = 1, 3 do
    angle = angles[index]
    open("gpu_age_" .. angle .. ".fem")
    mi_loadsolution()
    torque = mo_gapintegral("SlidingBand", 0)
    ia, va, lambda_a = mo_getcircuitproperties("A")
    ib, vb, lambda_b = mo_getcircuitproperties("B")
    ic, vc, lambda_c = mo_getcircuitproperties("C")
    gap_a_0 = mo_getgapa("SlidingBand", 0)
    h1_acc, h1_acs, h1_brc, h1_brs, h1_btc, h1_bts =
        mo_getgapharmonics("SlidingBand", 1)
    h6_acc, h6_acs, h6_brc, h6_brs, h6_btc, h6_bts =
        mo_getgapharmonics("SlidingBand", 6)
    print("XFEMM_GPU_AGE_RESULT," .. angle .. "," .. torque .. "," ..
          lambda_a .. "," .. lambda_b .. "," .. lambda_c .. "," .. gap_a_0 ..
          "," .. h1_brc .. "," .. h1_brs .. "," .. h1_btc .. "," .. h1_bts ..
          "," .. h6_brc .. "," .. h6_brs .. "," .. h6_btc .. "," .. h6_bts)
    mo_close()
    mi_close()
end

quit()
