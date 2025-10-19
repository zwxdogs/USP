// 线性阵列，继承自探头基类probe

#include "probe.H"

class linear_probe : public probe
{
private:
    double kerf;
    double pitch;
    int N_el;
    double width;   // x direction
    double height;  // y direction

public:
    // Constructor
    linear_probe(long fc, double kerf, double pitch, int N_el, double height, double bandwidth, double focus_z, double focus_x, double focus_y) :
        probe(fc, bandwidth, focus_z, focus_x, focus_y)
    {
        this->kerf = kerf;
        this->pitch = pitch;
        this->N_el = N_el;
        this->width = pitch - kerf;
        this->height = height;
    };

    // Destructor
    ~linear_probe() { };
};