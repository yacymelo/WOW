#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <cstdlib>
#include <stdlib.h>
#include <vector>

//#include <sys/time.h>
//#include <parallel/algorithm>
//#include <omp.h>
#include <cassert>

#include "Edgelist.h"
#include "Vault.h"
#include "Config.h"

using namespace std;


template <class V, class E=int>
class Graph {
    
public:
    int num_vertices;
    unsigned long num_edges;
    vector<Vault<V>> Vault_imp;
    V* vertices;
    V* vertices_tmp;
    unsigned long row_ratio[XBAR_SIZE];
    //bool vertexpropertyowner;
    //int tiles_per_dim;
    //int num_threads;
        
public:
    Graph(): num_vertices(0), num_edges(0), vertices(nullptr), vertices_tmp(nullptr) {
        for (auto i=0; i<NUM_VAULT; i++) {
            Vault<V> V_exa;
            Vault_imp.push_back(V_exa);
        }
        for (auto i=0; i<XBAR_SIZE; i++)
            row_ratio[i] = 0;
    }
    void LoadGraph(const char* filename);
    void GetRatio(vector<Vault<V>>* Vault_imp, unsigned long* row_ratio);
    void ParttoVault(edgelist_t<E>* A_edges, vector<Vault<V>>* Vault_imp, int Vault_num, int Xbar_size);
    void VertexInti(vector<Vault<V>>* Vault_imp, int Vault_num);
    void VertexInti_G();
    void SetActive(int v_id);
    void SetRoot(int v_id);
    bool HaveActive();
    void VetexIDtoVault(int v_id, int* vault_num, int* id_num);
    void SetActive_inVault(int v_id);
    void GRun_Syn_Pro(int iterations);
    void Vertex_Tmp_Inti_G();
    void ResetVertex();
    void ResetVertexTmp();

//    void ReadEdgelist(GraphMat::edgelist_t<E> A_edges);
//    void getVertexEdgelist(GraphMat::edgelist_t<V> & myedges);
//    void getEdgelist(GraphMat::edgelist_t<E> & myedges);
//    void ReadMTX(const char* filename);
//    void ReadGraphMatBin(const char* filename);
//    void WriteGraphMatBin(const char* filename);
//
//    void setAllActive();
//    void setAllInactive();
//    void setActive(int v);
//    void setInactive(int v);
//
//    void setAllVertexproperty(const V& val);
//    void setVertexproperty(int v, const V& val);
//    V getVertexproperty(int v) const;
//    bool vertexNodeOwner(const int v) const;
//    void saveVertexproperty(std::string fname, bool includeHeader=true) const;
//    void reset();
//    void shareVertexProperty(Graph<V,E>& g);
//    int getNumberOfVertices() const;
//    void applyToAllVertices(void (*ApplyFn)(const V&, V*, void*), void* param=nullptr);
//    template<class T> void applyReduceAllVertices(T* val, void (*ApplyFn)(V*, T*, void*), void (*ReduceFn)(const T&, const T&,T*,void*)=AddFn<T>, void* param=nullptr);
//    void applyToAllEdges(void (*ApplyFn)(E*, const V&, const V&, void*), void* param=nullptr);
//    ~Graph();
//
//private:
//    int vertexToNative(int vertex, int nsegments, int len) const;
//    int nativeToVertex(int vertex, int nsegments, int len) const;
    
};

template<class V, class E>
void Graph<V,E>::VertexInti(vector<Vault<V>>* Vault_imp, int Vault_num) {
    V vp;
    for (auto i=0; i<Vault_num; i++) {
        int star, v_num;
        star = (*Vault_imp)[i].star_ver;
        v_num = (*Vault_imp)[i].num_vertices;
        (*Vault_imp)[i].vertices = reinterpret_cast<V *>(malloc((size_t)v_num * sizeof(V)));
        for (auto id_v = 0; id_v<v_num; id_v++){
            (*Vault_imp)[i].vertices[id_v] = vp;
            (*Vault_imp)[i].vertices[id_v].id = id_v + star;
        }
    }
}

template<class V, class E>
void Graph<V,E>::VertexInti_G() {
    V vp;
    vertices = reinterpret_cast<V *>(malloc((size_t)num_vertices * sizeof(V)));
    for (auto i=0; i<num_vertices; i++){
        vertices[i] = vp;
        vertices[i].id = i;
    }
}

template<class V, class E>
void Graph<V,E>::Vertex_Tmp_Inti_G() {
    vertices_tmp = reinterpret_cast<V *>(malloc((size_t)num_vertices * sizeof(V)));
    for (auto i=0; i<num_vertices; i++)
        vertices_tmp[i] = vertices[i];
}

template<class V, class E>
void Graph<V,E>::ParttoVault(edgelist_t<E>* A_edges, vector<Vault<V>>* Vault_imp, int Xbar_size, int Vault_num) {
    int xbar_part;
    if (A_edges->n % Xbar_size) //ceil(xbar_part)
        xbar_part = A_edges->n / Xbar_size + 1;
    else
        xbar_part = A_edges->n / Xbar_size;
    unsigned long * xbar_part_num = reinterpret_cast<unsigned long *>(malloc((size_t)xbar_part * sizeof(unsigned long)));
    for (auto i=0; i<xbar_part; i++) //inti
        xbar_part_num[i] = 0;
    for (auto i=0; i<A_edges->num_edges; i++)
        xbar_part_num[A_edges->edges[i].dst / Xbar_size] ++;
    vector<edge_t<E> *> xbar_part_edges;
    for (auto i=0; i<xbar_part; i++) {
        xbar_part_edges.push_back(reinterpret_cast<edge_t<E> *>(malloc((uint64_t)xbar_part_num[i] * sizeof(edge_t<E>))));
    }
    for (auto i=0; i<xbar_part; i++) //inti
        xbar_part_num[i] = 0;
    for (auto i=0; i<A_edges->num_edges; i++){
        int j = A_edges->edges[i].dst / Xbar_size;
        auto m = xbar_part_num[j];
        xbar_part_edges[j][m].src = A_edges->edges[i].src;
        xbar_part_edges[j][m].dst = A_edges->edges[i].dst;
        xbar_part_edges[j][m].val = A_edges->edges[i].val;
        xbar_part_num[j] ++;
    }
    
//    for (auto i=0; i<xbar_part; i++) {
//        for(auto j=0; j<xbar_part_num[i]; j++){
//            cout << '(' << xbar_part_edges[i][j].src << ',' << xbar_part_edges[i][j].dst << ',' << xbar_part_edges[i][j].val << ')' << endl;
//        }
//        cout << "----------------------------" << endl;
//    }
    unsigned long vault_part_num[Vault_num];
    for (auto i=0; i<Vault_num; i++)
        vault_part_num[i] = 0;
    int v_x_times = xbar_part / Vault_num;
    int v_x_off = xbar_part % Vault_num;
    int index_v_x_star = 0;
    int index_v_x_end = 0;
    for (auto i=0; i<Vault_num; i++){
        index_v_x_star = index_v_x_end;
        for (auto j=0; j<v_x_times; j++){
            vault_part_num[i] = vault_part_num[i] + xbar_part_num[index_v_x_end];
            index_v_x_end ++;
        }
        if (v_x_off>i){
            vault_part_num[i] = vault_part_num[i] + xbar_part_num[index_v_x_end];
            index_v_x_end ++;
        }
        int cur_vault = i;
        
        (*Vault_imp)[cur_vault].edges = reinterpret_cast<edge_t<E> *>(malloc((size_t)vault_part_num[cur_vault] * sizeof(edge_t<E>)));
        
        unsigned long index_vault = 0;
        for (auto num_i=0; num_i<v_x_times; num_i++){
            for (auto j=0; j<xbar_part_num[index_v_x_star+num_i]; j++){
                (*Vault_imp)[cur_vault].edges[index_vault].src = xbar_part_edges[index_v_x_star+num_i][j].src;
                (*Vault_imp)[cur_vault].edges[index_vault].dst = xbar_part_edges[index_v_x_star+num_i][j].dst;
                (*Vault_imp)[cur_vault].edges[index_vault].val = xbar_part_edges[index_v_x_star+num_i][j].val;
                index_vault ++;
            }
            if (xbar_part_edges[index_v_x_star+num_i] != nullptr)
                delete xbar_part_edges[index_v_x_star+num_i];
            //delete xbar_part_edges[index_v_x_star+num_i];
        }
        if (v_x_off>cur_vault){
            for (auto j=0; j<xbar_part_num[index_v_x_star+v_x_times]; j++){
                (*Vault_imp)[cur_vault].edges[index_vault].src = xbar_part_edges[index_v_x_end-1][j].src;
                (*Vault_imp)[cur_vault].edges[index_vault].dst = xbar_part_edges[index_v_x_end-1][j].dst;
                (*Vault_imp)[cur_vault].edges[index_vault].val = xbar_part_edges[index_v_x_end-1][j].val;
                index_vault ++;
            }
            if (xbar_part_edges[index_v_x_end-1] != nullptr)
                delete xbar_part_edges[index_v_x_end-1];
        }
        (*Vault_imp)[cur_vault].num_edges = vault_part_num[cur_vault];
        (*Vault_imp)[cur_vault].star_ver = index_v_x_star * Xbar_size;
        (*Vault_imp)[cur_vault].end_ver = index_v_x_end * Xbar_size -1;
        (*Vault_imp)[cur_vault].num_vertices = (index_v_x_end - index_v_x_star)*Xbar_size;
    }
    
//    for (auto i=0; i<xbar_part; i++) {
//        delete xbar_part_edges[i];
//    }
    xbar_part_edges.clear();
    vector<edge_t<E> *>().swap(xbar_part_edges);
    free(xbar_part_num);
    xbar_part_num =nullptr;
    //for (auto i=0; i<Vault_num; i++){
}

template<class V, class E>
void Graph<V,E>::GetRatio(vector<Vault<V>>* Vault_imp, unsigned long* row_ratio){
    for(auto ind_va=0; ind_va<NUM_VAULT; ind_va++){
        unsigned long star = 0;
        int cur_src = (*Vault_imp)[ind_va].edges[star].src;
        int i;
        for(i=0; i<(*Vault_imp)[ind_va].num_edges; i++){
            //cout << cur_src;
            if(cur_src != (*Vault_imp)[ind_va].edges[i].src || (i-star) == XBAR_SIZE){
                row_ratio[i-star-1] ++;
                star = i;
                cur_src = (*Vault_imp)[ind_va].edges[star].src;
            }
        }
        row_ratio[i-star-1] ++;
    }
//    unsigned long star = 0;
//    int cur_src = A_edges->edges[star].src;
//    for(auto i=0; i<num_edges; i++){
//        cout << cur_src;
//        if(cur_src != A_edges->edges[i].src || i-star == XBAR_SIZE){
//            row_ratio[i-star] ++;
//            star = i;
//            cur_src = A_edges->edges[star].src;
//        }
//    }
}



template<class V, class E>
void Graph<V,E>::LoadGraph(const char* filename) {
    edgelist_t<E> A_edges;
    load_edgelist(filename, &A_edges);// binary format with header and edge weights
    num_vertices = A_edges.m + 1;
    num_edges = A_edges.num_edges;
    int Vault_num = NUM_VAULT;
    int Xbar_size = XBAR_SIZE;
    cout << "num_edges:" << num_edges << endl;
    ParttoVault(&A_edges, &Vault_imp, Xbar_size, Vault_num);
    GetRatio(&Vault_imp, row_ratio);
    unsigned long sub_numbers = 0;
    
    for(auto i=0; i<NUM_VAULT; i++) {
        unsigned long num_s = 0;
        //cout << "<<<<<" << endl;
        //int star_ver  = Vault_imp[i].star_ver;
        //int end_ver = Vault_imp[i].end_ver;
        for (auto x=Vault_imp[i].star_ver; x<Vault_imp[i].end_ver; ){
            for (auto y=0; y<num_vertices; ){

                //cout << "coming" << endl;
                bool index = true;
                bool index_2 = true;
                while (num_s < Vault_imp[i].num_edges &&  index_2) {
                    bool in_sub = false;
                    
                    if (Vault_imp[i].edges[num_s].src >= y && Vault_imp[i].edges[num_s].src < y+XBAR_SIZE)
                        if(Vault_imp[i].edges[num_s].dst >= x && Vault_imp[i].edges[num_s].dst < x+XBAR_SIZE)
                            in_sub = true;
                    //cout << Vault_imp[i].edges[num_s].src << " " << Vault_imp[i].edges[num_s].dst << endl;
                    if(in_sub){
                        if(index){
                            sub_numbers++;
                            index = false;
                            //cout << "____" << endl;
                        }
                        num_s ++;
                    }
                    else {
                        index_2 = false;
                    }
                }
                y= y + XBAR_SIZE;
            }
            x = x + XBAR_SIZE;
        }
        //cout << sub_numbers << endl;
    }
    
    cout <<sub_numbers << endl;
    cout << "=========" << endl;
    
    unsigned long sub_edges = 0;
    for(auto i=0; i<XBAR_SIZE; i++){
        sub_edges += row_ratio[i] * (i+1);
        cout << i+1 << ":" << row_ratio[i] << endl;
    }
    cout << "sub_edges:" << sub_edges << endl;
    //VertexInti(&Vault_imp, Vault_num);
    VertexInti_G();
        

//    for(auto i=0; i<Vault_num; i++){
//        cout << "star_ver:" << Vault_imp[i].star_ver << endl;
//        cout << "end_ver:" << Vault_imp[i].end_ver << endl;
//        cout << "num_vertices:" << Vault_imp[i].num_vertices << endl;
//        cout << "num_edges:" << Vault_imp[i].num_edges << endl;
//        for(auto j=0; j<Vault_imp[i].num_edges; j++){
//            cout << '(' << Vault_imp[i].edges[j].src << ',' << Vault_imp[i].edges[j].dst << ',' << Vault_imp[i].edges[j].val << ')' << endl;
//        }
//        cout << "------------------------------" << endl;
//    }
}

template<class V, class E>
void Graph<V,E>::VetexIDtoVault(int v_id, int* vault_num, int* id_num){
    for (auto i=0; i<NUM_VAULT; i++){
        if(v_id >= Vault_imp[i].star_ver && v_id <= Vault_imp[i].end_ver){
            *vault_num = i;
            *id_num = v_id - Vault_imp[i].star_ver;
            break;
        }
    }
}

template<class V, class E>
void Graph<V,E>::SetActive_inVault(int v_id) {
    int vault_num, id_num;
    VetexIDtoVault(v_id, &vault_num, &id_num);
    Vault_imp[vault_num].vertices[id_num].active = true;
}

template<class V, class E>
void Graph<V,E>::SetActive(int v_id) {
    vertices[v_id].active = true;
}

template<class V, class E>
void Graph<V,E>::SetRoot(int v_id) {
    vertices[v_id].RootInti();
}

template<class V, class E>
bool Graph<V,E>::HaveActive() {
//    for (auto i=0; i<NUM_VAULT; i++) {
//        for (auto j=0; j<Vault_imp[i].num_vertices; j++)
//            if(Vault_imp[i].vertices[j].active)
//                return true;
//    }
//    return false;
    for (auto i=0; i<num_vertices; i++) {
        if(vertices[i].active)
            return true;
    }
    return false;
}

template<class V, class E>
void Graph<V,E>::ResetVertex() {
    for (auto i=0; i<num_vertices; i++)
        vertices[i].Reset();
}

template<class V, class E>
void Graph<V,E>::ResetVertexTmp() {
    for (auto i=0; i<num_vertices; i++)
        vertices_tmp[i].ResetTmp();
}



template<class V, class E>
void Graph<V,E>::GRun_Syn_Pro(int iterations) {
    assert(iterations == -1 || iterations>0);
    Vertex_Tmp_Inti_G();
    int iter_nums = 0;
    if (iterations == UNTIL_CONVERGENCE) {
        while (HaveActive()) {
//            cout << "iter_nums:" << iter_nums << endl;
//            for (auto i=0; i<num_vertices; i++)
//                if (vertices[i].active)
//                    cout << vertices[i].id << '\t';
//            cout << endl;
            for (auto i=0; i<NUM_VAULT; i++)
                Vault_imp[i].VRun_Syn_Pro(vertices, vertices_tmp);
            ResetVertex();
            for (auto i=0; i<num_vertices; i++)
                vertices[i].Apply(vertices_tmp[i]);
            ResetVertexTmp();
            iter_nums ++;
        }
    } else if (iterations>0) {
        for (auto iter=0; iter<iterations; iter++) {
            for (auto i=0; i<NUM_VAULT; i++) {
                Vault_imp[i].VRun_Syn_Pro(vertices, vertices_tmp);
            }
            ResetVertex();
            for (auto i=0; i<num_vertices; i++)
                vertices[i].Apply(vertices_tmp[i]);
            ResetVertexTmp();
        }
    }
    cout << "iter_nums:" << iter_nums << endl;
    
//    for (auto i=0; i<num_vertices; i++){
//        cout << vertices[i].id << ':' << vertices[i].depth <<endl;
//    }
    unsigned long acc_sum = 0;
    for (auto i=0; i<NUM_VAULT; i++) {
        cout << i <<"_remote_ac:" << Vault_imp[i].remote_ac << endl;
        acc_sum += Vault_imp[i].remote_ac;
    }
    cout << "acc_sum:" << acc_sum << endl;
}
