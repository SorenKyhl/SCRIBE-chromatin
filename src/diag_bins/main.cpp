#include <iostream>
#include <vector>
#include <fstream>


std::vector<int> read_bins()
{
    std::vector<int> bins;
    std::ifstream file;
    file.open("bins.txt");

    int bin;
    while(file >> bin)
    {
        bins.push_back(bin);
    }

    file.close();
    return bins;

}


std::vector<int> set_diag_bins(std::vector<int> diag_bin_boundaries, int nbeads)
{
    std::vector<int> diag_bin_lookup = std::vector<int>(nbeads, 0);

    int curr = 0;
    int bin_id = 0;
    for (int i = 0; i < nbeads; i++)
    {
        if (i >= diag_bin_boundaries[curr])
        {
            bin_id += 1;
            curr += 1;
        }
        std::cout << i << " " << bin_id << std::endl;
        diag_bin_lookup[i] = bin_id;
    }
    return diag_bin_lookup;
}


int main()
{
    std::vector<int> bins = read_bins();
    std::cout << bins[2] << std::endl;

    int nbeads = 10;
    std::vector<int> diag_bin_lookup = set_diag_bins(bins, nbeads);
    
    for (int e : diag_bin_lookup)
    {
         std::cout << e << " ";
    }
}
