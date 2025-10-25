# Calculate spatial impulse response.


def calc_SIR(probe, scan, para):
    # Get required parameters
    probe_line = probe.get_line()
    probe_corners = probe.get_corners()
    scan_position = scan.get_scan_xyz()
    dt = para.get_dt()
    c0 = para.get_c0()

    pass
    # To be implemented
