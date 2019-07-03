using namespace std;

template <typename T>
struct edge_t {
    edge_t() {}
    edge_t(int _src, int _dst, T _val)
    {
        src = _src;
        dst = _dst;
        val = _val;
    }
    int src;
    int dst;
    T val;
};

template <typename T>
struct edgelist_t {
    edge_t<T>* edges;
    int m;
    int n;
    unsigned long num_edges;
    edgelist_t() : m(0), n(0), num_edges(0), edges(nullptr) {}
    edgelist_t(int _m, int _n, int _num_edges)
    {
        m = _m;
        n = _n;
        num_edges = _num_edges;
        if(num_edges > 0) {
            edges = reinterpret_cast<edge_t<T>*>(malloc((size_t)num_edges * sizeof(edge_t<T>)));
        }
    }
    edgelist_t(edge_t<T>* edges, int m, int n, int num_edges) : edges(edges), m(m), n(n), num_edges(num_edges) {}
    void clear() {
        if (num_edges > 0) {
            _mm_free(edges);
        }
        edges = nullptr;
        num_edges = 0;
        m = 0;
        n = 0;
    }
};

template<typename T>
bool readLine (FILE * ifile, int * src, int * dst, T * val, bool binaryformat=false, bool edgeweights=false)
{
    if(binaryformat) {
        auto fread_bytes = fread(src, sizeof(int), 1, ifile);
        if (feof(ifile)) return false;
        assert(fread_bytes == 1);
        fread_bytes = fread(dst, sizeof(int), 1, ifile);
        if (feof(ifile)) return false;
        assert(fread_bytes == 1);
        if (edgeweights) {
            fread_bytes = fread(val, sizeof(T), 1, ifile);
            if (feof(ifile)) return false;
            assert(fread_bytes == 1);
        } else {
            *val = (T)(1);
        }
    } else {
        if (edgeweights) {
            int ret;
            if (std::is_same<T, float>::value) {
                ret = fscanf(ifile, "%d %d %f", src, dst, val);
                if (ret != 3) return false;
            } else if (std::is_same<T, double>::value) {
                ret = fscanf(ifile, "%d %d %lf", src, dst, val);
                if (ret != 3) return false;
            } else if (std::is_same<T, int>::value) {
                ret = fscanf(ifile, "%d %d %d", src, dst, val);
                if (ret != 3) return false;
            } else if (std::is_same<T, unsigned int>::value) {
                ret = fscanf(ifile, "%d %d %u", src, dst, val);
                if (ret != 3) return false;
            }else {
                std::cout << "Data type not supported (read)" << std::endl;
            }
        } else {
            int ret = fscanf(ifile, "%d %d", src, dst);
            if (ret == 2) {
                *val = (T)(1);
            } else return false;
        }
        if (feof(ifile)) return false;
    }
    return true;
}

template<typename T>
void get_maxid_and_num_edges(FILE* fp, int* m, int* n, unsigned long int* num_edges, bool binaryformat=false, bool header=false, bool edgeweights=false) {
    if (header) {
        int tmp_[3];
        if (binaryformat) {
            auto fread_bytes = fread(tmp_, sizeof(int), 3, fp);
            assert(fread_bytes == 3);
            *m = tmp_[0];
            *n = tmp_[1];
            *num_edges = tmp_[2];
        } else {
            int ret = fscanf(fp, "%d %d %u", &(tmp_[0]), &(tmp_[1]), &(tmp_[2]));
            assert(ret == 3);
            *m = tmp_[0];
            *n = tmp_[1];
            *num_edges = tmp_[2];
        }
        return;
    } else { //no header
        unsigned long num_edges_ = 0;
        int tempsrc, tempdst;
        int maxm = 0;
        int maxn = 0;
        T tempval;
        while(true) {
            if(feof(fp)) {
                break;
            }
            if (!readLine<T>(fp, &tempsrc, &tempdst, &tempval, binaryformat, edgeweights)) {
                break;
            }
            maxm = (maxm > tempsrc)?(maxm):(tempsrc);
            maxn = (maxn > tempdst)?(maxn):(tempdst);
            num_edges_++;
        }
        *m = maxm;
        *n = maxn;
        *num_edges = num_edges_;
    }
}

template <typename T>
void load_edgelist(const char* dir, edgelist_t<T>* edgelist, bool binaryformat=false, bool header=false, bool edgeweights=false) {
    edgelist->m = 0;
    edgelist->n = 0;
    edgelist->num_edges = 0;
    std::stringstream fname_ss;
    fname_ss << dir;

    FILE* fp;
    if (binaryformat) {
        fp = fopen(fname_ss.str().c_str(), "rb");
    } else {
        fp = fopen(fname_ss.str().c_str(), "r");
    }
    if(!fp) {
        printf("Could not open file: %s\n", fname_ss.str().c_str());
        exit(-1);
    } else {
        printf("Reading file: %s\n", fname_ss.str().c_str());
    }
    int m_, n_;
    unsigned long num_edges_;
    get_maxid_and_num_edges<T>(fp, &m_, &n_, &num_edges_, binaryformat, header, edgeweights);
    int vertex_max_;
    vertex_max_ = std::max(m_, n_);
    edgelist->m = std::max(vertex_max_, edgelist->m);
    edgelist->n = std::max(vertex_max_, edgelist->n);
    edgelist->num_edges = num_edges_;
    fclose(fp);
    
    edgelist->edges = reinterpret_cast<edge_t<T>*>(malloc((uint64_t)edgelist->num_edges * (uint64_t)sizeof(edge_t<T>)));
    
    if (binaryformat) {
        fp = fopen(fname_ss.str().c_str(), "rb");
    } else {
        fp = fopen(fname_ss.str().c_str(), "r");
    }
    if(!fp) {
        exit(-1);
    }
    if (header) { //remove header
        int m_, n_;
        unsigned long num_edges_;
        get_maxid_and_num_edges<T>(fp, &m_, &n_, &num_edges_, binaryformat, header, edgeweights);
    }
    
    unsigned long int cnt = 0;
    while(true) {
        if (feof(fp)) {
            break;
        }
        if (!readLine<T>(fp, &(edgelist->edges[cnt].src), &(edgelist->edges[cnt].dst), &(edgelist->edges[cnt].val), binaryformat, edgeweights)) {
            break;
        }
#ifdef __DEBUG
        //std::cout <<(edgelist->edges[nnzcnt].src) << " " << (edgelist->edges[nnzcnt].dst) << std::endl;
        if(edgelist->edges[cnt].src <= 0 ||
           edgelist->edges[cnt].dst <= 0 ||
           edgelist->edges[cnt].src > edgelist->m ||
           edgelist->edges[cnt].dst > edgelist->n)
        {
            std::cout << "Invalid edge, i, j, nnz: " << i << " , " << j << " , " << nnzcnt << std::endl;
            exit(0);
        }
#endif
        cnt++;
    }
    fclose(fp);
//    for (auto i = 0; i < cnt; i++)
//        cout << '(' << edgelist->edges[i].src << ',' << edgelist->edges[i].dst << ',' << edgelist->edges[i].val << ')' << endl;
}
