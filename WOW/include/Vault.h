#include "Config.h"

using namespace std;

template <class V, class E=int>
class Vault {
    
public:
    
    int star_ver;
    int end_ver;
    int num_vertices;
    unsigned long num_edges;
    V* vertices;
    V* vertices_tmp;
    edge_t<E>* edges;
    unsigned long remote_ac = 0;
    //bool vertexpropertyowner;
    //int tiles_per_dim;
    //int num_threads;
    
public:
    Vault(): star_ver(0), end_ver(0), num_vertices(0), num_edges(0), vertices(nullptr), vertices_tmp(nullptr), edges(nullptr){}
    void Inti_tmp();
    void GetXbarLine(int* g_src, int* g_num_edges, unsigned long Runed_edges);
    void VRun_Syn_Pro(V* vertices_G, V* vertices_G_tmp);
    
};

template<class V, class E>
void Vault<V,E>::Inti_tmp() {
    vertices_tmp = reinterpret_cast<V *>(malloc((size_t)num_vertices * sizeof(V)));
    for (auto i=0; i<num_vertices; i++)
        vertices_tmp[i] = vertices[i];
}

template<class V, class E>
void Vault<V,E>::GetXbarLine(int* g_src, int* g_num_edges, unsigned long Runed_edges) {
    *g_src = edges[Runed_edges].src;
    int i;
    for (i=0; i<(num_edges - Runed_edges) && i<XBAR_SIZE; i++){
        if (*g_src != edges[i+Runed_edges].src){
            *g_num_edges = i;
            return;
        }
    }
    if (i == (num_edges - Runed_edges) || i == XBAR_SIZE){
        *g_num_edges = i;
    }
}

template<class V, class E>
void Vault<V,E>::VRun_Syn_Pro(V* vertices_G, V* vertices_G_tmp) {
    unsigned long Runed_edges = 0;
    while (Runed_edges < num_edges) {
        int g_src, g_num_edges;
        GetXbarLine(&g_src, &g_num_edges, Runed_edges);
//        cout << "active:" << vertices_G[14].active << endl;
//        cout << "g_src:" << g_src << endl;
//        cout << "g_num_edges:" << g_num_edges << endl;
//        cout << "active:" << vertices_G[g_src].active << endl;
        if(vertices_G[g_src].active) {
            if(g_src < star_ver || g_src > end_ver)
                remote_ac ++;
            for (auto i_edge=0; i_edge<g_num_edges; i_edge++){
                V Res;
                int Res_dst = edges[Runed_edges+i_edge].dst;
                vertices_G[g_src].Process_Edge(Res);
                vertices_G_tmp[Res_dst].Reduce(Res);
            }
        }
        Runed_edges = Runed_edges + g_num_edges;
    }
}

