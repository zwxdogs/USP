# Calculate spatial impulse response.
import simulation.tools.SIR_calc as SIR_cpp


def calc_SIR(probe, scan, para):
    # Call C++ extension for calculation
    match probe.get_symbol():
        case "rect":
            SIR_result = SIR_cpp.calc_polygon(
                probe.get_line(),
                probe.get_calc_line(),
                probe.get_corners(),
                probe.get_calc_corners(),
                scan.get_scan_xyz(),
                para.get_dt(),
                para.get_c0(),
            )

        case _:
            raise ValueError("Unsupported probe symbol for SIR calculation.")

    return SIR_result
