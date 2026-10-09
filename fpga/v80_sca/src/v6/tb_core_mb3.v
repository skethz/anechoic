// r3 testbench (copy of tb_core_mb.v with the BIAS and HBX parameters of sca_core_mb_r3.v).
// Stream testbench for sca_core_mb.v (copy of tb_core.v): drives stim.hex (random input gaps), checks every output word
// against expect.hex (random output back-pressure). Prints the cycle of each trial's start, result and last word so per-trial
// cycles can be compared with v6.4b. Added for MB: larger stimulus memory (K coupling planes), a run-time check of the
// residue property the single-copy banks rely on (extractor e only issues rows x with x mod 8 == e), and a count of
// issued flips per extractor.
`timescale 1ns/1ps
module tb;
    parameter integer DRAIN = 2;
    parameter integer K = 2;
    parameter integer MAXCYC = 40000000;
    parameter integer BIAS = 0;
    parameter integer HBX = 0;
    reg clk = 1'b0, rst_n = 1'b0;
    always #1 clk = ~clk;

    reg [127:0] stim [0:524287];
    reg [63:0] expv [0:131071];
    integer nbeats, nwords, fd, rc;

    wire [127:0] s_tdata; wire s_tvalid, s_tready;
    wire [63:0] m_tdata; wire m_tvalid; reg m_tready;
    reg gap;
    integer bi, ei, errors, cycle, seedv, res_err, ee;

    sca_core_mb #(.K(K), .DRAIN(DRAIN), .BIAS(BIAS), .HBX(HBX)) dut (.clk(clk), .rst_n(rst_n), .s_axis_tdata(s_tdata), .s_axis_tvalid(s_tvalid),
        .s_axis_tready(s_tready), .m_axis_tdata(m_tdata), .m_axis_tvalid(m_tvalid), .m_axis_tready(m_tready));

    assign s_tdata = stim[bi];
    assign s_tvalid = (bi < nbeats) & ~gap;

    initial begin
        $readmemh("stim.hex", stim);
        $readmemh("expect.hex", expv);
        fd = $fopen("counts.txt", "r");
        rc = $fscanf(fd, "%d %d", nbeats, nwords);
        $fclose(fd);
        $display("beats %0d words %0d DRAIN %0d K %0d", nbeats, nwords, DRAIN, K);
        bi = 0; ei = 0; errors = 0; cycle = 0; seedv = 12345; gap = 1'b0; m_tready = 1'b1; res_err = 0;
        repeat (20) @(posedge clk);
        @(negedge clk) rst_n = 1'b1;   // release between edges (no race with the DUT's synchronous reset)
    end

    always @(posedge clk) begin
        cycle <= cycle + 1;
        if (rst_n) begin
            if (s_tvalid & s_tready) bi <= bi + 1;
            gap <= ($random(seedv) % 16) == 0;
            m_tready <= ($random(seedv) % 8) != 0;
            if (m_tvalid & m_tready) begin
                if (m_tdata !== expv[ei]) begin
                    errors = errors + 1;
                    if (errors <= 20) $display("MISMATCH word %0d: got %016h expected %016h (cycle %0d)", ei, m_tdata, expv[ei], cycle);
                end
                ei <= ei + 1;
                if (ei + 1 == nwords) begin
                    $display("DONE: %0d words, %0d errors, %0d residue violations, cycle %0d", nwords, errors, res_err, cycle);
                    // (xsim prints a string-valued conditional as a number, so two separate $display calls)
                    if (errors == 0 && res_err == 0) $display("TESTBENCH PASSED"); else $display("TESTBENCH FAILED");
                    $finish;
                end
            end
        end
        if (cycle == MAXCYC) begin $display("TIMEOUT at word %0d beat %0d", ei, bi); $display("TESTBENCH FAILED"); $finish; end
    end

    // residue property: a valid flip issued by extractor e has x mod 8 == e (the bank of port e holds only those rows)
    always @(posedge clk) if (rst_n)
        for (ee = 0; ee < 8; ee = ee + 1)
            if (dut.av_bus[ee] && (dut.ax_bus[11 * ee +: 3] != ee)) begin
                res_err = res_err + 1;
                if (res_err <= 5) $display("RESIDUE VIOLATION extractor %0d x %0d (cycle %0d)", ee, dut.ax_bus[11 * ee +: 11], cycle);
            end

    // trial timing: entry to SEED (trial start) and to OUT (results ready), from the DUT state register
    reg [3:0] st_prev; initial st_prev = 4'hf;
    always @(posedge clk) begin
        st_prev <= dut.st;
        if (dut.st == 4'd5 && st_prev != 4'd5) $display("trial start: cycle %0d", cycle);
        if (dut.st == 4'd9 && st_prev != 4'd9) $display("trial out: cycle %0d", cycle);
    end
    // launch timing: header accepted -> last word of the launch accepted
    always @(posedge clk) if (rst_n & s_tvalid & s_tready & (dut.st == 4'd0))
        $display("launch start: beat %0d cycle %0d", bi, cycle);
    // per-trial timing: report when the 34th word (flips) of a trial is accepted
    integer wi_in_trial; initial wi_in_trial = 0;
    always @(posedge clk) if (rst_n & m_tvalid & m_tready) begin
        if (wi_in_trial == 33) $display("trial end: word %0d cycle %0d flips %0d", ei, cycle, m_tdata);
        wi_in_trial <= wi_in_trial + 1;
    end
    // trial word counting depends on trace length; reset the counter on the testbench side from the DUT's FSM
    always @(posedge clk) if (dut.st == 4'd4) wi_in_trial <= 0;
endmodule
