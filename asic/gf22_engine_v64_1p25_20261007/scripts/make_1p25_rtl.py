#!/usr/bin/env python3
"""Derive the 1.25 GHz variant from the v6.4b ASIC RTL (1 GHz memory variant): the step-table SRAM read becomes a
two-cycle path. The table macro output is captured two cycles after its read edge (enable tab_ld3) and the next-step
registers load one cycle later than before (tab_ld4 / tab_ld5). The table entry for the next step is read at the start
of the current step and consumed only after that step's decision batch (>= ~17 cycles later), and reads are >= 8 cycles
apart, so no cycle count changes (verified against the unmodified v6.4b RTL). Synthesis gets a matching multicycle
constraint (setup 2 / hold 1) from the table macro to tab_q."""
import sys

src, dst = sys.argv[1], sys.argv[2]
s = open(src).read()


def rep(old, new):
    global s
    assert s.count(old) == 1, old[:70]
    s = s.replace(old, new)


rep("""// sca_core_asic.v -- ASIC variant (GF 22FDX estimate) of the V80 Onsager-SCA engine core v6.4b (sca_core.v, 8 extractors).""",
    """// sca_core_asic_1p25.v -- 1.25 GHz ASIC variant: as sca_core_asic.v (v6.4b memories), plus a two-cycle step-table read
// (the table macro output is captured two cycles after its read edge; see "1.25 GHz" notes). Cycle counts unchanged.
// sca_core_asic.v -- ASIC variant (GF 22FDX estimate) of the V80 Onsager-SCA engine core v6.4b (sca_core.v, 8 extractors).""")
rep("""    reg tab_ld3, tab_ld2, tab_ld1, tab_ld4;    // table read in flight: address with tab_ld1, data in tab_q with tab_ld3""",
    """    reg tab_ld3, tab_ld2, tab_ld1, tab_ld4, tab_ld5;  // table read in flight: address with tab_ld1, data in tab_q with tab_ld4""")
rep("""    always @(posedge clk) tab_q <= tab_q1;""",
    """    // 1.25 GHz: capture the macro output two cycles after the read edge (multicycle path 2: the output holds until the
    // next read, which is >= 8 cycles later)
    always @(posedge clk) if (tab_ld3) tab_q <= tab_q1;""")
rep("""        tab_ld2 <= tab_ld1; tab_ld3 <= tab_ld2; tab_ld4 <= tab_ld3;
        if (tab_ld3) begin
            tn_fourT <= tab_q[31:0]; tn_q <= tab_q[63:32]; tn_kconst <= tab_q[95:64]; tn_kcorr <= tab_q[127:96];
        end
        if (tab_ld4) begin tn_U <= tn_twoT - 32'd1; tn_L <= 32'd1 - tn_twoT; end""",
    """        tab_ld2 <= tab_ld1; tab_ld3 <= tab_ld2; tab_ld4 <= tab_ld3; tab_ld5 <= tab_ld4;
        if (tab_ld4) begin   // 1.25 GHz: one cycle later than v6.4b (still >= 13 cycles before first use)
            tn_fourT <= tab_q[31:0]; tn_q <= tab_q[63:32]; tn_kconst <= tab_q[95:64]; tn_kcorr <= tab_q[127:96];
        end
        if (tab_ld5) begin tn_U <= tn_twoT - 32'd1; tn_L <= 32'd1 - tn_twoT; end""")
open(dst, 'w').write(s)
print('written', dst)
# Note: the generated file's default CPL_G was then set to 2 (two groups per 256 x 128 coupling macro) by a second
# asserted replacement (see the CPL_G parameter comment in rtl/sca_core_asic_1p25.v).
