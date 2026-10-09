#!/usr/bin/env python3
"""Package the verified sca_v80 RTL with a passive start-to-done cycle counter (same mechanism as v80_snowball).

Timer registers are placed at free AXI-Lite addresses above the generated argument map. Their offsets are written to
results/instrumentation.json and to build/sca_timer_offsets.h, which the host compiles against.
"""
from pathlib import Path
import hashlib, json, os, re, zipfile
import xml.etree.ElementTree as ET

root = Path(os.environ['TASK_ROOT'])
variant = os.environ.get('SCA_VARIANT', '')
sfx = f'_{variant}' if variant else ''
proj = f'hls_sca{sfx}'
assert (root / f'results/hls{sfx}_exit_status').read_text().strip() == '0', 'HLS job must finish successfully'
assert (root / f'results/cosim{sfx}_exit_status').read_text().strip() == '0', 'RTL co-simulation job must pass'
assert '*** C/RTL co-simulation finished: PASS ***' in (root / f'logs/cosim{sfx}.log').read_text(), 'RTL verification is required'
archive = root / f'build/{proj}/solution/impl/ip/xilinx_com_hls_sca_v80_1_0.zip'
iprepo = Path(os.environ.get('SNOWBALL_IP_REPO', str(root / 'build/iprepo')))
assert iprepo.resolve().is_relative_to(root / 'build'), 'IP repository must stay in the scratch build directory'
dest = iprepo / 'sca_v80'
dest.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(archive) as z:
    for name in z.namelist():
        p = Path(name)
        if p.is_absolute() or '..' in p.parts:
            raise ValueError(name)
    z.extractall(dest)
p = dest / 'hdl/verilog/sca_v80_control_s_axi.v'
text = p.read_text()
original_hash = hashlib.sha256(p.read_bytes()).hexdigest()

# Generated map: lines like "ADDR_S_DATA_0 = 8'h70,". Choose three free word addresses above the highest used one.
addrs = [int(v, 16) for v in re.findall(r"ADDR_\w+\s*=\s*\d+'h([0-9a-fA-F]+)", text)]
width = int(re.search(r"ADDR_BITS\s*=\s*(\d+)", text).group(1))
top = max(addrs)
ident, lo, hi = top + 0x10, top + 0x14, top + 0x18
if hi >= (1 << width):
    # Small register maps (v5.2: 7 address bits): use the lowest unused word addresses of the generated map instead.
    # Small register maps (v5.2: 7 address bits, every word named): reuse per-argument ADDR_*_CTRL placeholder words that
    # the generated RTL defines but never decodes (the parameter name occurs only at its definition).
    ctrl = [(name, int(v, 16)) for name, v in re.findall(r"(ADDR_\w+_CTRL)\s*=\s*\d+'h([0-9a-fA-F]+)", text)
            if name != 'ADDR_AP_CTRL' and len(re.findall(rf"\b{name}\b", text)) == 1]
    assert len(ctrl) >= 3, f'no three undecoded placeholder words below 2^{width}'
    ident, lo, hi = sorted(a for _, a in ctrl)[:3]
    addrs = [a for a in addrs if a not in (ident, lo, hi)]
assert hi < (1 << width), f'no free register space below 2^{width}'
for a in (ident, lo, hi):
    assert f"{width}'h{a:x}:" not in text.lower() and a not in addrs

counter = '''
// Passive timing monitor (copied from v80_snowball): first asserted ap_start to ap_done, host auto-restart disabled.
reg perf_active = 1'b0;
reg [63:0] perf_counter = 64'd0;
reg [63:0] perf_completed = 64'd0;
always @(posedge ACLK) begin
    if (ARESET) begin
        perf_active <= 1'b0;
        perf_counter <= 64'd0;
        perf_completed <= 64'd0;
    end else begin
        if (!perf_active && ACLK_EN && ap_start) begin
            perf_active <= 1'b1;
            perf_counter <= 64'd0;
        end else if (perf_active) begin
            perf_counter <= perf_counter + 64'd1;
            if (ap_done) begin
                perf_completed <= perf_counter + 64'd1;
                perf_active <= 1'b0;
            end
        end
    end
end
'''
marker = '//------------------------Instantiation------------------'
assert text.count(marker) == 1
text = text.replace(marker, counter + '\n' + marker)
marker = '            case (raddr)'
assert text.count(marker) == 1
text = text.replace(marker, marker + f'''
                {width}'h{ident:x}: rdata <= 32'h53434131; // SCA1: timing monitor present
                {width}'h{lo:x}: rdata <= perf_completed[31:0];
                {width}'h{hi:x}: rdata <= perf_completed[63:32];''')
p.write_text(text)

component = dest / 'component.xml'
for _, (prefix, uri) in ET.iterparse(component, events=['start-ns']):
    ET.register_namespace(prefix, uri)
tree = ET.parse(component)
ns = {'s': 'http://www.spiritconsortium.org/XMLSchema/SPIRIT/1685-2009'}
for path, tag in [('s:model/s:views', 's:view'), ('s:fileSets', 's:fileSet')]:
    parent = tree.getroot().find(path, ns)
    for child in list(parent):
        if (child.findtext('s:name', namespaces=ns) or '').startswith('xilinx_vhdl'):
            parent.remove(child)
tree.write(component, encoding='UTF-8', xml_declaration=True)

(root / f'build/sca_timer_offsets{sfx}.h').write_text(
    f'#pragma once\n#define SCA_TIMER_ID 0x{ident:x}\n#define SCA_CYC_LO 0x{lo:x}\n#define SCA_CYC_HI 0x{hi:x}\n')
(root / f'results/instrumentation{sfx}.json').write_text(json.dumps({
    'ip_archive_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
    'control_rtl_original_sha256': original_hash,
    'control_rtl_instrumented_sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
    'addr_bits': width, 'identity_offset': ident, 'cycles_low_offset': lo, 'cycles_high_offset': hi,
    'identity_value': '0x53434131', 'counter_bits': 64,
    'counter_scope': 'Kernel-clock edges from first ap_start to ap_done; includes table loading, (optional) matrix staging, '
                     'all trials of the launch and output writes; excludes host and PCIe staging.',
    'kernel_datapath_modified_by_instrumentation': False,
}, indent=2) + '\n')
print(dest, hex(ident), hex(lo), hex(hi))
