// class probe for all of the probe.

class probe
{
private:
    long fc;
    double kerf;
    double pitch;
    int N_el;
    double width;
    double height;
    double bandwidth;
    double focus_z;
    double focus_x;
    double focus_y;
    double center[3];

public:
    // Constructor
    probe(long fc, double kerf, double pitch, int N_el, double height, double bandwidth, double focus_z, double focus_x, double focus_y)
    {
        this->fc = fc;
        this->kerf = kerf;
        this->pitch = pitch;
        this->N_el = N_el;
        this->width = pitch - kerf;
        this->height = height;
        this->bandwidth = bandwidth;
        this->focus_z = focus_z;
        this->focus_x = focus_x;
        this->focus_y = focus_y;
        this->center[0] = 0;
        this->center[1] = 0;
        this->center[2] = 0;
    };

    // Destructor
    ~probe();
};