#include <armadillo>
#include <cmath>
#include <numeric>

using namespace arma;
using namespace std;

//===============================================
// Basis generation using the lexicographic order
//===============================================

void generate_basis(int M, int N, Mat<int>& basis) {
    vector<int> state(M, 0);
    state[0] = N;
    
    // Initialize basis matrix with proper dimensions
    basis.set_size(0, M);

    while (true) {
        basis.insert_rows(basis.n_rows, 1);
        for (int i = 0; i < M; i++) basis(basis.n_rows - 1, i) = state[i];

        int k = M - 2;
        while (k >= 0 && state[k] == 0) --k;
        if (k < 0) break;

        state[k]--;
        state[k+1] = N - accumulate(state.begin(), state.begin() + k + 1, 0);
        for (int i = k+2; i < M; i++) state[i] = 0;
    }
}

//===============================================
// Unique Tag 
//===============================================

double tag_vector(const Row<int>& vec, const vector<int>& primes) {
    double tag = 0.0;
    for (size_t i = 0; i < vec.n_elem; i++)
        tag += sqrt(primes[i]) * vec[i];
    return tag;
}

//===============================================
// Binary Search 
//===============================================

uword find_index_by_tag(double tag, const vec& tags, const uvec& sort_index, double eps = 1e-12) {
    // manual binary search over the *sorted* view (tags(sort_index))
    uword lo = 0, hi = sort_index.n_elem;
    while (lo < hi) {
        uword mid = (lo + hi) / 2;
        double tmid = tags(sort_index(mid));
        if (tmid < tag) lo = mid + 1;
        else hi = mid;
    }
    if (lo >= sort_index.n_elem) throw runtime_error("Tag not found (upper bound).");
    double found = tags(sort_index(lo));
    if (std::abs(found - tag) > eps) throw runtime_error("Tag not found (tolerance).");
    return sort_index(lo);
}

//===============================================
// 1D Bose Hubbard Hamiltonian
//===============================================

sp_mat build_hamiltonian(const Mat<int>& basis,
    int M, int N, double J, double U,
    const vector<int>& primes,
    const vec& tags, const uvec& sort_index)
{
const uword D = basis.n_rows;

std::vector<uword> I; I.reserve(D * (M * 2 + 1));
std::vector<uword> Jc; Jc.reserve(D * (M * 2 + 1));
std::vector<double> V; V.reserve(D * (M * 2 + 1));

auto right = [&](int i) { return (i + 1) % M; }; 

for (uword v = 0; v < D; ++v) {
const Row<int> state = basis.row(v);

// Diagonal
double diag = 0.0;
for (int i = 0; i < M; ++i) {
int n = state(i);
diag += 0.5 * U * n * (n - 1);
}
I.push_back(v); Jc.push_back(v); V.push_back(diag);

// Kinetic
for (int i = 0; i < M; ++i) {
int j = right(i);

// Hopping from j to i 
if (state(j) > 0) {
Row<int> s2 = state;
s2(i) += 1;
s2(j) -= 1;
double amp = -J * std::sqrt( static_cast<double>(state(i) + 1) * static_cast<double>(state(j)) );
double t = tag_vector(s2, primes);
uword u = find_index_by_tag(t, tags, sort_index);
I.push_back(u); Jc.push_back(v); V.push_back(amp);
}

// Hopping from i to j 
if (state(i) > 0) {
Row<int> s2 = state;
s2(j) += 1;
s2(i) -= 1;
double amp = -J * std::sqrt( static_cast<double>(state(j) + 1) * static_cast<double>(state(i)) );
double t = tag_vector(s2, primes);
uword u = find_index_by_tag(t, tags, sort_index);
I.push_back(u); Jc.push_back(v); V.push_back(amp);
}
}
}

// Build sparse matrix 
// Create locations matrix where each column is [row, col]
umat locations(2, I.size());
vec  vals(V.size());
for (uword k = 0; k < I.size(); ++k) {
locations(0, k) = I[k];   // row index
locations(1, k) = Jc[k]; // col index
vals(k) = V[k];
}
sp_mat H(locations, vals, D, D, /*sort_locations=*/true, /*check_for_zeros=*/false);
return H;
}


int main() {


int M = 4, N = 4;
double J = 1.0, U = 2.0;

Mat<int> basis;
generate_basis(M, N, basis);
cout << "Basis size D = " << basis.n_rows << endl;


vector<int> primes = {3, 5, 7, 11};  
vec tags(basis.n_rows);
for (uword v = 0; v < basis.n_rows; v++)
tags(v) = tag_vector(basis.row(v), primes);
uvec sort_idx = sort_index(tags);     


sp_mat H = build_hamiltonian(basis, M, N, J, U, primes, tags, sort_idx);


vec evals;
mat evecs;
eigs_sym(evals, evecs, H, /*k=*/2, "sa");

evals.print("Lowest eigenvalues (E0, E1):");

vec gs = evecs.col(0);

return 0;
}
