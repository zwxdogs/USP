// 所有探头的基类probe

class probe
{
private:
    long fc;
    double bandwidth;
    double focus_z;
    double focus_x;
    double focus_y;
    double center[3];

public:
    // Constructor
    probe(long fc, double bandwidth, double focus_z, double focus_x, double focus_y)
    {
        this->fc = fc;
        this->bandwidth = bandwidth;
        this->focus_z = focus_z;
        this->focus_x = focus_x;
        this->focus_y = focus_y;
        this->center[0] = 0;
        this->center[1] = 0;
        this->center[2] = 0;
    };

    // Destructor
    ~probe() { };
};