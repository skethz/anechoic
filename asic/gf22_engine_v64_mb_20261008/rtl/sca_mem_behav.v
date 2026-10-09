// Behavioral models of the single-port SRAM macros used by sca_core_mb_asic.v (multi-bit ASIC variant; as the v6.4b ASIC,
// plus the bias memory). Same ports and timing
// as the synthesis wrappers around the compiler macros (single-port, flow-through, high-density bitcell):
//   all inputs sampled at the rising clock edge; ce = 1 selects an access, we = 1 a write (all bits written);
//   read (ce & ~we): rdata = mem[addr] after the edge (the macro output acts as the read register);
//   write (ce & we): mem[addr] = wdata, rdata keeps its previous value; NOP (~ce): rdata keeps its previous value.
// Synthesis replaces mem_cpl / mem_tab / mem_trc by wrappers that instantiate the compiler macros; those wrappers are
// kept on the synthesis server only (the macro views are under NDA).
`timescale 1ns/1ps
`default_nettype none

module mem_sp_behav #(
    parameter integer WORDS = 256,
    parameter integer W = 64,
    parameter integer AW = 8
) (
    input  wire          clk,
    input  wire          ce,
    input  wire          we,
    input  wire [AW-1:0] addr,
    input  wire [W-1:0]  wdata,
    output reg  [W-1:0]  rdata
);
    reg [W-1:0] mem [0:WORDS-1];
    always @(posedge clk) begin
        if (ce) begin
            if (we) mem[addr] <= wdata;
            else rdata <= mem[addr];
        end
    end
endmodule

// coupling residue-bank slice: 256 rows x W bits (W = 64 * CPL_G)
module mem_cpl #(
    parameter integer W = 256
) (
    input  wire         clk,
    input  wire         ce,
    input  wire         we,
    input  wire [7:0]   addr,
    input  wire [W-1:0] wdata,
    output wire [W-1:0] rdata
);
    mem_sp_behav #(.WORDS(256), .W(W), .AW(8)) m (.clk(clk), .ce(ce), .we(we), .addr(addr), .wdata(wdata), .rdata(rdata));
endmodule

// per-step table: 4096 x 128
module mem_tab (
    input  wire         clk,
    input  wire         ce,
    input  wire         we,
    input  wire [11:0]  addr,
    input  wire [127:0] wdata,
    output wire [127:0] rdata
);
    mem_sp_behav #(.WORDS(4096), .W(128), .AW(12)) m (.clk(clk), .ce(ce), .we(we), .addr(addr), .wdata(wdata), .rdata(rdata));
endmodule

// per-step n_lin trace: 4096 x 12
module mem_trc (
    input  wire         clk,
    input  wire         ce,
    input  wire         we,
    input  wire [11:0]  addr,
    input  wire [11:0]  wdata,
    output wire [11:0]  rdata
);
    mem_sp_behav #(.WORDS(4096), .W(12), .AW(12)) m (.clk(clk), .ce(ce), .we(we), .addr(addr), .wdata(wdata), .rdata(rdata));
endmodule
// bias memory: 256 x 128 (2048 x int16); same macro type as a CPL_G = 2 coupling slice
module mem_bias (
    input  wire         clk,
    input  wire         ce,
    input  wire         we,
    input  wire [7:0]   addr,
    input  wire [127:0] wdata,
    output wire [127:0] rdata
);
    mem_cpl #(.W(128)) m (.clk(clk), .ce(ce), .we(we), .addr(addr), .wdata(wdata), .rdata(rdata));
endmodule
`default_nettype wire
