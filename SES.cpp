#include "SES.h"

using namespace std;


// This routine finds a set of belt patch with two (or one) neighboring concave patches, all has areas below 
// the min_area as defined below. then it merge those patches and mark it as a concave patch. all the neighbor 
// lists are updated as required.
void 
MySES_Raw :: melt_patch_belt()
{
  double  tmp_min_area_belt = 1;
  double  tmp_min_area_concave = 2.0;
  //double    tmp_min_area_convex = 0.5;

  size_t i, j;

  vector<MySES_Patch*>                tmp_patch_toDel;
  vector<MySES_Patch*>                tmp_patch_toAdd;

  tmp_patch_toDel.clear();    
  tmp_patch_toAdd.clear();


  vector<MySES_Patch*>::iterator      it_patch;

  for (i = 0; i < m_patch_all.size(); i ++){

    if (m_patch_all[i]->m_patch_type != 1) continue;     // Not a belt type patch, so ignore
    if (m_patch_all[i]->m_area > tmp_min_area_belt) continue; // area greater than min-belt area, so ignore

    //now the current patch is a belt and it has an area greater than the min
    if (m_patch_all[i]->m_nb_concave.size() == 2) {     // usually a belt patch has two concave neighboring patch.
      MySES_Patch*            tmp_concave_first = m_patch_all[i]->m_nb_concave[0];   // 1st convace neighbor patch, concave1
      MySES_Patch*            tmp_concave_second = m_patch_all[i]->m_nb_concave[1];  // 2nd concave neighbor patch, concave2
      MySES_Patch*            tmp_belt_current = m_patch_all[i];  // this variable holds the current patch, belt

      //both neighbor concave patches are small, merge them
      if (tmp_concave_first->m_area < tmp_min_area_concave &&
        tmp_concave_second->m_area < tmp_min_area_concave) {
        set<MySES_Patch*>   tmp_external_all;
        tmp_external_all.clear();

        //deleting the 1st concave as the belts neighbor
        if (tmp_belt_current->m_nb_concave.size() > 0){
          it_patch = find(tmp_belt_current->m_nb_concave.begin(), tmp_belt_current->m_nb_concave.end(), tmp_concave_first);
          if (it_patch != tmp_belt_current->m_nb_concave.end()){
            tmp_belt_current->m_nb_concave.erase(it_patch); // erasing concave1 from concave neighbor list of belt
          }
        }

        //deleting the 2nd concave as the belts neighbor
        if (tmp_belt_current->m_nb_concave.size() > 0){
          it_patch = find(tmp_belt_current->m_nb_concave.begin(), tmp_belt_current->m_nb_concave.end(), tmp_concave_second);
          if (it_patch != tmp_belt_current->m_nb_concave.end()){
            tmp_belt_current->m_nb_concave.erase(it_patch); // erasing concave2 from concave neighbor list of belt
          }
        }

        // erasing belt from belt neighbor list of the concave1
        if (tmp_concave_first->m_nb_belt.size() > 0){
          it_patch = find(tmp_concave_first->m_nb_belt.begin(), tmp_concave_first->m_nb_belt.end(), tmp_belt_current);
          if (it_patch != tmp_concave_first->m_nb_belt.end()){
            tmp_concave_first->m_nb_belt.erase(it_patch);
          }
        }

        // erasing concave2 from concave neighbor list of the concave1    
        if (tmp_concave_first->m_nb_concave.size() > 0){
          it_patch = find(tmp_concave_first->m_nb_concave.begin(), tmp_concave_first->m_nb_concave.end(), tmp_concave_second);
          if (it_patch != tmp_concave_first->m_nb_concave.end()){
            tmp_concave_first->m_nb_concave.erase(it_patch);
          }
        }

        // erasing belt from belt neighbor list of the concave2
        if (tmp_concave_second->m_nb_belt.size() > 0){
          it_patch = find(tmp_concave_second->m_nb_belt.begin(), tmp_concave_second->m_nb_belt.end(), tmp_belt_current);
          if (it_patch != tmp_concave_second->m_nb_belt.end()){
            tmp_concave_second->m_nb_belt.erase(it_patch);
          }
        }

        // erasing concave1 from concave neighbor list of the concave2       
        if (tmp_concave_second->m_nb_concave.size() > 0){
          it_patch = find(tmp_concave_second->m_nb_concave.begin(), tmp_concave_second->m_nb_concave.end(), tmp_concave_first);
          if (it_patch != tmp_concave_second->m_nb_concave.end()){
            tmp_concave_second->m_nb_concave.erase(it_patch);
          }
        }

        //create the new patch of type cocave
        MySES_Patch*    tmp_new_concave = new MySES_Patch();
        tmp_new_concave->m_patch_type = 2;  // Making a new patch of concave type
        tmp_new_concave->m_area = tmp_concave_first->m_area + tmp_concave_second->m_area + tmp_belt_current->m_area;
        tmp_new_concave->m_verts.clear();  // Vertex list is empty now

        // some sanity check
        if (tmp_concave_first->m_verts.size() == 0){
          cout << "ERROR! tmp_concave_first->m_verts.size() == 0 !" << endl;
        }

        if (tmp_concave_second->m_verts.size() == 0){
          cout << "ERROR! tmp_concave_second->m_verts.size() == 0 !" << endl;
        }

        if (tmp_belt_current->m_verts.size() == 0){
          cout << "ERROR! tmp_belt_current->m_verts.size() == 0 !" << endl;
        }

        if (tmp_concave_first->m_triangles.size() == 0){
          cout << "ERROR! tmp_concave_first->m_triangles.size() == 0 !" << endl;
        }

        if (tmp_concave_second->m_triangles.size() == 0){
          cout << "ERROR! tmp_concave_second->m_triangles.size() == 0 !" << endl;
        }

        if (tmp_belt_current->m_triangles.size() == 0){
          cout << "ERROR! tmp_belt_current->m_triangles.size() == 0 !" << endl;
        }

        // inserting all the vertices of these three patches, 1 belt and 2 concaves
        tmp_new_concave->m_verts.insert(tmp_concave_first->m_verts.begin(), tmp_concave_first->m_verts.end());
        tmp_new_concave->m_verts.insert(tmp_concave_second->m_verts.begin(), tmp_concave_second->m_verts.end());
        tmp_new_concave->m_verts.insert(tmp_belt_current->m_verts.begin(), tmp_belt_current->m_verts.end());

        // inserting all the triangles of these three patches, 1 belt and 2 concaves
        tmp_new_concave->m_triangles.insert(tmp_new_concave->m_triangles.end(), tmp_concave_first->m_triangles.begin(), tmp_concave_first->m_triangles.end());
        tmp_new_concave->m_triangles.insert(tmp_new_concave->m_triangles.end(), tmp_concave_second->m_triangles.begin(), tmp_concave_second->m_triangles.end());
        tmp_new_concave->m_triangles.insert(tmp_new_concave->m_triangles.end(), tmp_belt_current->m_triangles.begin(), tmp_belt_current->m_triangles.end());

        //collect the external patches (basically all the neighboring patches of the above 3 patches.
        tmp_external_all.insert(tmp_belt_current->m_nb_belt.begin(), tmp_belt_current->m_nb_belt.end());
        tmp_external_all.insert(tmp_belt_current->m_nb_concave.begin(), tmp_belt_current->m_nb_concave.end());
        tmp_external_all.insert(tmp_belt_current->m_nb_convex.begin(), tmp_belt_current->m_nb_convex.end());
        tmp_external_all.insert(tmp_concave_first->m_nb_belt.begin(), tmp_concave_first->m_nb_belt.end());
        tmp_external_all.insert(tmp_concave_first->m_nb_concave.begin(), tmp_concave_first->m_nb_concave.end());
        tmp_external_all.insert(tmp_concave_first->m_nb_convex.begin(), tmp_concave_first->m_nb_convex.end());
        tmp_external_all.insert(tmp_concave_second->m_nb_belt.begin(), tmp_concave_second->m_nb_belt.end());
        tmp_external_all.insert(tmp_concave_second->m_nb_concave.begin(), tmp_concave_second->m_nb_concave.end());
        tmp_external_all.insert(tmp_concave_second->m_nb_convex.begin(), tmp_concave_second->m_nb_convex.end());


        set<MySES_Patch*>::iterator it_ext;
        //delete the external links
        for (it_ext = tmp_external_all.begin(); it_ext != tmp_external_all.end(); it_ext ++) { //for each neighbor patch
          if ((*it_ext)->m_nb_belt.size() > 0){  
            it_patch = find((*it_ext)->m_nb_belt.begin(), (*it_ext)->m_nb_belt.end(), tmp_belt_current);
            if (it_patch != (*it_ext)->m_nb_belt.end()){
              (*it_ext)->m_nb_belt.erase(it_patch); // delete the current belt from the belt-list
            }
          }
          if ((*it_ext)->m_nb_concave.size() > 0){
            it_patch = find((*it_ext)->m_nb_concave.begin(), (*it_ext)->m_nb_concave.end(), tmp_concave_first);
            if (it_patch != (*it_ext)->m_nb_concave.end()){
              (*it_ext)->m_nb_concave.erase(it_patch); // delete the concave1 from the concave-list
            }
          }
          if ((*it_ext)->m_nb_concave.size() > 0){
            it_patch = find((*it_ext)->m_nb_concave.begin(), (*it_ext)->m_nb_concave.end(), tmp_concave_second);
            if (it_patch != (*it_ext)->m_nb_concave.end()){
              (*it_ext)->m_nb_concave.erase(it_patch); // delete the concave2 from the concave-list
            }
          }
        }

        //setup the new links
        for (it_ext = tmp_external_all.begin(); it_ext != tmp_external_all.end(); it_ext ++){ //for each neighboring patch
          (*it_ext)->m_nb_concave.push_back(tmp_new_concave);  // newly created patch goes into their concave-neighbor list
          if ((*it_ext)->m_patch_type == 1){
            tmp_new_concave->m_nb_belt.push_back((*it_ext)); //newly created patch gets ins belt-neighbor list filled
          }
          else if ((*it_ext)->m_patch_type == 2){
            tmp_new_concave->m_nb_concave.push_back((*it_ext)); //newly created patch gets ins concave-neighbor list filled
          }
          else if ((*it_ext)->m_patch_type == 3){
            tmp_new_concave->m_nb_convex.push_back((*it_ext)); //newly created patch gets ins convex-neighbor list filled
          }
        }

        tmp_patch_toAdd.push_back(tmp_new_concave);   // adding the new patch to list of patch to add
        tmp_patch_toDel.push_back(tmp_concave_first); // adding the old 3 patches to the list of patch to delete
        tmp_patch_toDel.push_back(tmp_concave_second);
        tmp_patch_toDel.push_back(tmp_belt_current);
      }

      else {
        MySES_Patch* tmp_small_concave = NULL;
        MySES_Patch* tmp_big_concave = NULL;

         //only the first neighbor concave patch is small
        if (tmp_concave_first->m_area < tmp_min_area_concave){
          tmp_small_concave = tmp_concave_first;
          tmp_big_concave = tmp_concave_second;
        }

        //only the second neighbor concave patch is small
        else if (tmp_concave_second->m_area < tmp_min_area_concave){
          tmp_small_concave = tmp_concave_second;
          tmp_big_concave = tmp_concave_first;
        }

        //both are big, skip;
        if (tmp_small_concave == NULL || tmp_big_concave == NULL) continue;

        //break the internal links
        if (tmp_belt_current->m_nb_concave.size() > 0){
          it_patch = find(tmp_belt_current->m_nb_concave.begin(), tmp_belt_current->m_nb_concave.end(), tmp_small_concave);
          if (it_patch != tmp_belt_current->m_nb_concave.end()){
            tmp_belt_current->m_nb_concave.erase(it_patch); // removing the convace_small from the belt's neighbor list
          }
        }

        if (tmp_small_concave->m_nb_belt.size() > 0){
          it_patch = find(tmp_small_concave->m_nb_belt.begin(), tmp_small_concave->m_nb_belt.end(), tmp_belt_current);
          if (it_patch != tmp_small_concave->m_nb_belt.end()){
            tmp_small_concave->m_nb_belt.erase(it_patch); // removing the belt from the concave_small's neighbor list
          }
        }

        //create the new patch, patch type is concave
        MySES_Patch*    tmp_new_patch = new MySES_Patch();
        tmp_new_patch->m_patch_type = 2;
        tmp_new_patch->m_area = tmp_small_concave->m_area + tmp_belt_current->m_area;

        // some sanity check
        if (tmp_small_concave->m_triangles.size() == 0){
          cout << "ERROR: tmp_small_concave->m_triangles.size() == 0" << endl;
        }

        if (tmp_belt_current->m_triangles.size() == 0){
          cout << "ERROR: tmp_belt_current->m_triangles.size() == 0" << endl;
        }

        if (tmp_small_concave->m_verts.size() == 0){
          cout << "ERROR: tmp_small_concave->m_verts.size() == 0" << endl;
        }

        if (tmp_belt_current->m_verts.size() == 0){
          cout << "ERROR: tmp_belt_current->m_verts.size() == 0" << endl;
        }

        tmp_new_patch->m_triangles.insert(tmp_new_patch->m_triangles.end(), tmp_small_concave->m_triangles.begin(), tmp_small_concave->m_triangles.end());
        tmp_new_patch->m_triangles.insert(tmp_new_patch->m_triangles.end(), tmp_belt_current->m_triangles.begin(), tmp_belt_current->m_triangles.end());
        tmp_new_patch->m_verts.insert(tmp_small_concave->m_verts.begin(), tmp_small_concave->m_verts.end());
        tmp_new_patch->m_verts.insert(tmp_belt_current->m_verts.begin(), tmp_belt_current->m_verts.end());

        //collect the external patches
        set<MySES_Patch*>       tmp_external_all;
        tmp_external_all.clear();

        tmp_external_all.insert(tmp_small_concave->m_nb_convex.begin(), tmp_small_concave->m_nb_convex.end());
        tmp_external_all.insert(tmp_small_concave->m_nb_concave.begin(), tmp_small_concave->m_nb_concave.end());
        tmp_external_all.insert(tmp_small_concave->m_nb_belt.begin(), tmp_small_concave->m_nb_belt.end());
        tmp_external_all.insert(tmp_belt_current->m_nb_convex.begin(), tmp_belt_current->m_nb_convex.end());
        tmp_external_all.insert(tmp_belt_current->m_nb_concave.begin(), tmp_belt_current->m_nb_concave.end());
        tmp_external_all.insert(tmp_belt_current->m_nb_belt.begin(), tmp_belt_current->m_nb_belt.end());

        set<MySES_Patch*>::iterator it_ext;
        //break the external links
        for (it_ext = tmp_external_all.begin(); it_ext != tmp_external_all.end(); it_ext ++){
          if ((*it_ext)->m_nb_concave.size() > 0){
            it_patch = find((*it_ext)->m_nb_concave.begin(), (*it_ext)->m_nb_concave.end(), tmp_small_concave);
            if (it_patch != (*it_ext)->m_nb_concave.end()){
              (*it_ext)->m_nb_concave.erase(it_patch);
            }
          }

          if ((*it_ext)->m_nb_belt.size() > 0){
            it_patch = find((*it_ext)->m_nb_belt.begin(), (*it_ext)->m_nb_belt.end(), tmp_belt_current);
            if (it_patch != (*it_ext)->m_nb_belt.end()){
              (*it_ext)->m_nb_belt.erase(it_patch);
            }
          }
        }

        //setup the new links
        for (it_ext = tmp_external_all.begin(); it_ext != tmp_external_all.end(); it_ext ++){
          (*it_ext)->m_nb_concave.push_back(tmp_new_patch);
          if ((*it_ext)->m_patch_type == 1){
            tmp_new_patch->m_nb_belt.push_back((*it_ext));
          }
          else if ((*it_ext)->m_patch_type == 2){
            tmp_new_patch->m_nb_concave.push_back((*it_ext));
          }
          else if ((*it_ext)->m_patch_type == 3){
            tmp_new_patch->m_nb_convex.push_back((*it_ext));
          }
        }
        tmp_patch_toDel.push_back(tmp_belt_current);
        tmp_patch_toDel.push_back(tmp_small_concave);
        tmp_patch_toAdd.push_back(tmp_new_patch);
      }
    }

  }// end of for over m_patch_all

  set<MySES_Patch*>   set_A_all(m_patch_all.begin(), m_patch_all.end());
  set<MySES_Patch*>   set_B_toDel(tmp_patch_toDel.begin(), tmp_patch_toDel.end());
  set<MySES_Patch*>   set_C_remain;

  //cout << "set_A_all: " << set_A_all.size() << "  " << m_patch_all.size() << endl;
  //cout << "set_B_toDel: " << set_B_toDel.size() << "  " << tmp_patch_toDel.size() << endl;

  set_difference(set_A_all.begin(), set_A_all.end(), set_B_toDel.begin(), set_B_toDel.end(), inserter(set_C_remain, set_C_remain.begin()));

  //cout << "set_C_remain: " << set_C_remain.size() << endl;

  m_patch_all.clear();
  m_patch_all.insert(m_patch_all.end(), set_C_remain.begin(), set_C_remain.end());
  m_patch_all.insert(m_patch_all.end(), tmp_patch_toAdd.begin(), tmp_patch_toAdd.end());

  //cout << "tmp_patch_toAdd: " << tmp_patch_toAdd.size() << endl;
  //cout << "m_patch_all: " << m_patch_all.size() << endl;

  /* Zujun had the following delete commented, I think it should not be commented */
  for (j = 0; j < tmp_patch_toDel.size(); j ++){
    delete tmp_patch_toDel[j];
  }

}

// This routine creates MySES_Patch vectors that contains all the patches of this surface
// it also initializes m_nb_belt, m_nb_concave, m_nb_convex of these patches and initialize
// m_patch_all (vector of a pointers to all patches)
// Also compute the area of a triangular face.
void 
MySES_Raw :: melt_patch_init() {
  size_t i;
  MySES_Patch* tmp_patch;

  // Indexing patch by integer, by mapping int-->patch pointer
  map<unsigned int, MySES_Patch*>             tmp_patch_all;

  //map patch id to type id, this can find which patch is of which type;
  map<unsigned int, unsigned int>             tmp_patchid_typeid;

  // for all the triangular faces of this surface
  for (i = 0; i < m_face_patch_id.size(); i ++){
    tmp_patchid_typeid[m_face_patch_id[i]] = m_face_type[i];  // mapping face_id --> face_type
  }

  map<unsigned int, set<TID_VERT> >::iterator it_patch;

  //create a MySES_Patch for each patch
  //starting will all the belt patches
  for (it_patch = m_patch_belt.begin(); it_patch != m_patch_belt.end(); it_patch ++){
    tmp_patch = new MySES_Patch(); 
    tmp_patch->m_patch_type = 1;
    tmp_patch->m_verts = it_patch->second;
    tmp_patch_all.insert(pair<unsigned int, MySES_Patch*>(it_patch->first, tmp_patch)); // inserting belt patches in all patch list
  }

  //Now the concave patches
  for (it_patch = m_patch_concave.begin(); it_patch != m_patch_concave.end(); it_patch ++){
    tmp_patch = new MySES_Patch();
    tmp_patch->m_patch_type = 2;
    tmp_patch->m_verts = it_patch->second;
    tmp_patch_all.insert(pair<unsigned int, MySES_Patch*>(it_patch->first, tmp_patch)); //inserting concave pathes in all patch list
  }

  //Now the convex patches
  for (it_patch = m_patch_convex.begin(); it_patch != m_patch_convex.end(); it_patch ++){
    tmp_patch = new MySES_Patch();
    tmp_patch->m_patch_type = 3;
    tmp_patch->m_verts = it_patch->second;
    tmp_patch_all.insert(pair<unsigned int, MySES_Patch*>(it_patch->first, tmp_patch)); //inserting convex pathes in all patch list
  }

  map<pair<unsigned int, unsigned int>, set<unsigned int> > tmp_edgeid_patchid;// mapping from edge to set of patches the edge belong
  tmp_edgeid_patchid.clear();

  unsigned int    edge_v1, edge_v2;

  //find patch id for each edge of all the triangles of this surface
  for (i = 0; i < m_face_triangle.size(); i ++){
    // each edge are like, v1---v2, where v1 always has lower id
    // and v2 has the higher id
    // A triangle has 3 edges, so three blocks below, for each edge, we are inserting the patch-id in the corresponding 
    // patch-set, for which this edge belongs. Note that, more than one patch can be inserted for an edge, so the second
    // data-type of tmp_edgeid_patch_id is a set
    edge_v1 = (m_face_triangle[i][0] < m_face_triangle[i][1])?   m_face_triangle[i][0] : m_face_triangle[i][1];
    edge_v2 = (m_face_triangle[i][0] < m_face_triangle[i][1])?   m_face_triangle[i][1] : m_face_triangle[i][0];
    tmp_edgeid_patchid[pair<unsigned int, unsigned int>(edge_v1, edge_v2)].insert(m_face_patch_id[i]); //create the mapping from edge-->patchid

    edge_v1 = (m_face_triangle[i][1] < m_face_triangle[i][2])?   m_face_triangle[i][1] : m_face_triangle[i][2];
    edge_v2 = (m_face_triangle[i][1] < m_face_triangle[i][2])?   m_face_triangle[i][2] : m_face_triangle[i][1];
    tmp_edgeid_patchid[pair<unsigned int, unsigned int>(edge_v1, edge_v2)].insert(m_face_patch_id[i]); //create the mapping from edge-->patchid

    edge_v1 = (m_face_triangle[i][0] < m_face_triangle[i][2])?   m_face_triangle[i][0] : m_face_triangle[i][2];
    edge_v2 = (m_face_triangle[i][0] < m_face_triangle[i][2])?   m_face_triangle[i][2] : m_face_triangle[i][0];
    tmp_edgeid_patchid[pair<unsigned int, unsigned int>(edge_v1, edge_v2)].insert(m_face_patch_id[i]); //create the mapping from edge-->patchid

    tmp_patch_all[m_face_patch_id[i]]->m_triangles.push_back(m_face_triangle[i]); // the vertex ids of triangle i is iserted under triangle-list of 
                                                                                  // the corresponding patch.
  }

  // defining some iterator variables to walk over the vertex list, triangle list etc.
  map<unsigned int, MySES_Patch*>::iterator   it_spatch;  
  set<TID_VERT>::iterator                     it_vert;
  vector<vector<TID_VERT> >::iterator         it_triangle;
  set<unsigned int>::iterator                 it_pid;

  //find neighbor patches 
  for (it_spatch = tmp_patch_all.begin(); it_spatch != tmp_patch_all.end(); it_spatch ++){ // for each patch
    set<unsigned int>   tmp_nb_patchid;  // set of patch-ids
    tmp_nb_patchid.clear();
    // all patches that are in patch-list of any of the triangle-edges of this patch are neighbor of this patch
    for (i = 0; i < it_spatch->second->m_triangles.size(); i ++) {  // for all the triangles of this patch
      edge_v1 = (it_spatch->second->m_triangles[i][0] < it_spatch->second->m_triangles[i][1])?
         		it_spatch->second->m_triangles[i][0]:it_spatch->second->m_triangles[i][1];

      edge_v2 = (it_spatch->second->m_triangles[i][0] < it_spatch->second->m_triangles[i][1])? 
			it_spatch->second->m_triangles[i][1]:it_spatch->second->m_triangles[i][0];
      tmp_nb_patchid.insert(tmp_edgeid_patchid[pair<unsigned int, unsigned int>(edge_v1, edge_v2)].begin(), tmp_edgeid_patchid[pair<unsigned int, unsigned int>(edge_v1, edge_v2)].end());
                
      edge_v1 = (it_spatch->second->m_triangles[i][1] < it_spatch->second->m_triangles[i][2])?
			it_spatch->second->m_triangles[i][1]:it_spatch->second->m_triangles[i][2];

      edge_v2 = (it_spatch->second->m_triangles[i][1] < it_spatch->second->m_triangles[i][2])? 
			it_spatch->second->m_triangles[i][2]:it_spatch->second->m_triangles[i][1];
      tmp_nb_patchid.insert(tmp_edgeid_patchid[pair<unsigned int, unsigned int>(edge_v1, edge_v2)].begin(), tmp_edgeid_patchid[pair<unsigned int, unsigned int>(edge_v1, edge_v2)].end());

      edge_v1 = (it_spatch->second->m_triangles[i][2] < it_spatch->second->m_triangles[i][0])?
			it_spatch->second->m_triangles[i][2]:it_spatch->second->m_triangles[i][0];

      edge_v2 = (it_spatch->second->m_triangles[i][2] < it_spatch->second->m_triangles[i][0])? 
			it_spatch->second->m_triangles[i][0]:it_spatch->second->m_triangles[i][2];
      tmp_nb_patchid.insert(tmp_edgeid_patchid[pair<unsigned int, unsigned int>(edge_v1, edge_v2)].begin(), tmp_edgeid_patchid[pair<unsigned int, unsigned int>(edge_v1, edge_v2)].end());
    }

    tmp_nb_patchid.erase(it_spatch->first); // a patch is not neighbor to itself, so deleting it from it's neighbor list

    for (it_pid = tmp_nb_patchid.begin(); it_pid != tmp_nb_patchid.end(); it_pid ++){
      if (tmp_patchid_typeid[*it_pid] == 1){ // neighbor-patch is belt
        it_spatch->second->m_nb_belt.push_back(tmp_patch_all[*it_pid]);
      }
      else if (tmp_patchid_typeid[*it_pid] == 2){ // neighbor-patch is concave
        it_spatch->second->m_nb_concave.push_back(tmp_patch_all[*it_pid]);
      }
      else if (tmp_patchid_typeid[*it_pid] == 3){ // neighbor-patch is convex
        it_spatch->second->m_nb_convex.push_back(tmp_patch_all[*it_pid]);
      }
    }
  }

  // Now finding area of a patch, a patch is a triangle.
  double tmp_area = 0;
  TID_VERT    vid_1;
  TID_VERT    vid_2;
  TID_VERT    vid_3;

  m_patch_all.clear();
  //compute the area for each patch
  for (it_spatch = tmp_patch_all.begin(); it_spatch != tmp_patch_all.end(); it_spatch ++){
    tmp_area = 0;
    for (i = 0; i < it_spatch->second->m_triangles.size(); i ++) { // for each triangles of this patch
      vid_1 = it_spatch->second->m_triangles[i][0];  // id of the 1st vertex
      vid_2 = it_spatch->second->m_triangles[i][1];  // id of the 2nd vertex
      vid_3 = it_spatch->second->m_triangles[i][2];  // id of the 3rd vertex

      // m_vert_xyz[vid] returns the 3D co-ordinate of the vertex with id = vid
      tmp_area += compute_triangle_area(m_vert_xyz[vid_1], m_vert_xyz[vid_2], m_vert_xyz[vid_3]);
    }
    it_spatch->second->m_area = tmp_area;
    m_patch_all.push_back(it_spatch->second);
  }
}

// This routine writes the patch_info of a surface, in 3 files one for each type of surface
void 
MySES_Raw ::write_patch_vtk(string pm_filename) {

  //cout << "m_patch_belt.size() = " << m_patch_belt.size() << endl;
  //cout << "m_patch_concave.size() = " << m_patch_concave.size() << endl;
  //cout << "m_patch_convex.size() = " << m_patch_convex.size()  << endl;

  write_patch_vtk(pm_filename + "_patchBelt.vtk", 1);
  write_patch_vtk(pm_filename + "_patchConcave.vtk", 2);
  write_patch_vtk(pm_filename + "_patchConvex.vtk", 3);
}

// This routine writes the specific type patch_info of surface to a file
void 
MySES_Raw :: write_patch_vtk(string pm_filename, int pm_patch_type) {
  size_t m;

  set<TID_VERT>               tmp_points_nodup;
  vector<vector<TID_VERT> >   tmp_triangles;  // will hold all the triangles that belongs to the patches of this type
	                                      // inner vector for vertex-id of 3 vertices of a triangle

  tmp_points_nodup.clear();
  tmp_triangles.clear();

  for (m = 0; m < m_patch_all.size(); m ++){
    //cout << "m = " << m << " out of " << m_patch_all.size() << endl;
    if (m_patch_all[m]->m_patch_type != pm_patch_type) continue;  // this patch is not of the type as arg-2 of this method
    if (m_patch_all[m]->m_triangles.size() > 0)
      tmp_triangles.insert(tmp_triangles.end(), m_patch_all[m]->m_triangles.begin(), m_patch_all[m]->m_triangles.end());
    if (m_patch_all[m]->m_verts.size() > 0)
      tmp_points_nodup.insert(m_patch_all[m]->m_verts.begin(), m_patch_all[m]->m_verts.end());
  }


  stringstream out_buffer;
  out_buffer << "# vtk DataFile Version 3.0" << endl;
  out_buffer << pm_filename << endl;
  out_buffer << "ASCII" << endl;
  out_buffer << "DATASET POLYDATA" << endl;

  vector<vector<double> > tmp_points;

  set<TID_VERT>::iterator it_vert;
  map<TID_VERT, int>  tmp_vert2id;

  for (m = 0, it_vert = tmp_points_nodup.begin(); it_vert != tmp_points_nodup.end(); m ++, it_vert ++){
    tmp_vert2id[(*it_vert)] = m;
    tmp_points.push_back(m_vert_xyz[(*it_vert)]);
  }

  out_buffer << "POINTS " << tmp_points.size() << " float" << endl;

  for(m = 0; m < tmp_points.size(); m ++){
    out_buffer << tmp_points[m][0] << " "
               << tmp_points[m][1] << " "
               << tmp_points[m][2] << endl;
  }

  out_buffer << "POLYGONS " << tmp_triangles.size() << " " << tmp_triangles.size() * 4 << endl;
  for (size_t i = 0; i < tmp_triangles.size(); i++){
    out_buffer << "3 ";
    out_buffer << tmp_vert2id[tmp_triangles[i][0]] << " ";
    out_buffer << tmp_vert2id[tmp_triangles[i][1]] << " ";
    out_buffer << tmp_vert2id[tmp_triangles[i][2]] << endl;
  }

  ofstream out_file(pm_filename.c_str());
  if (!out_file.is_open()){
    cerr << "ERROR: cannot open file " << pm_filename << endl;
    exit(1);
  }

  out_file << out_buffer.str();
  out_file.close();

}

// write all patches in a single file
void 
MySES_Raw :: write_patch_vtk_single(MySES_Patch* pm_patch, string pm_filename) {
  size_t m;

  set<TID_VERT>               tmp_points_nodup;
  for (m = 0; m < pm_patch->m_triangles.size(); m ++){
    tmp_points_nodup.insert(pm_patch->m_triangles[m][0]);
    tmp_points_nodup.insert(pm_patch->m_triangles[m][1]);
    tmp_points_nodup.insert(pm_patch->m_triangles[m][2]);
  }

  stringstream out_buffer;
  out_buffer << "# vtk DataFile Version 3.0" << endl;
  out_buffer << pm_filename << endl;
  out_buffer << "ASCII" << endl;
  out_buffer << "DATASET POLYDATA" << endl;

  vector<vector<double> > tmp_points;
  set<TID_VERT>::iterator it_vert;
  map<TID_VERT, int>  tmp_vert2id;

  for (m = 0, it_vert = tmp_points_nodup.begin(); it_vert != tmp_points_nodup.end(); m ++, it_vert ++){
    tmp_vert2id[(*it_vert)] = m;
    tmp_points.push_back(m_vert_xyz[(*it_vert)]);
  }

  out_buffer << "POINTS " << tmp_points.size() << " float" << endl;

  for(m = 0; m < tmp_points.size(); m ++){
    out_buffer << tmp_points[m][0] << " "
               << tmp_points[m][1] << " "
                << tmp_points[m][2] << endl;
  }

  out_buffer << "POLYGONS " << pm_patch->m_triangles.size() << " " << pm_patch->m_triangles.size() * 4 << endl;
  for (size_t i = 0; i < pm_patch->m_triangles.size(); i++){
    out_buffer << "3 ";
    out_buffer << tmp_vert2id[pm_patch->m_triangles[i][0]] << " ";
    out_buffer << tmp_vert2id[pm_patch->m_triangles[i][1]] << " ";
    out_buffer << tmp_vert2id[pm_patch->m_triangles[i][2]] << endl;
  }

  ofstream out_file(pm_filename.c_str());
    if (!out_file.is_open()){
      cerr << "ERROR: cannot open file " << pm_filename << endl;
      exit(1);
    }

  out_file << out_buffer.str();
  out_file.close();

}

// this routine writes the vertices as 3 vtk files 
void 
MySES_Raw ::write_vertex_vtk(string pm_filename) {
  write_vertex_vtk(pm_filename + "_vertBelt.vtk", 1);
  write_vertex_vtk(pm_filename + "_vertConcave.vtk", 2);
  write_vertex_vtk(pm_filename + "_vertConvex.vtk", 3);
}

 
// this routine writes specific type of vertices as one vtk files 
void 
MySES_Raw :: write_vertex_vtk(string pm_filename, int pm_vert_type) {
  stringstream out_buffer;
  out_buffer << "# vtk DataFile Version 3.0" << endl;
  out_buffer << pm_filename << endl;
  out_buffer << "ASCII" << endl;
  out_buffer << "DATASET POLYDATA" << endl;

  vector<vector<double> > tmp_points;
  size_t m;
  for (m = 0; m < m_vert_xyz.size(); m ++){
    if (m_vert_type[m] == pm_vert_type){
      tmp_points.push_back(m_vert_xyz[m]);
    }
  }

  out_buffer << "POINTS " << tmp_points.size() << " float" << endl;

  for(m = 0; m < tmp_points.size(); m ++){
    out_buffer << tmp_points[m][0] << " "
                << tmp_points[m][1] << " "
                << tmp_points[m][2] << endl;
  }

  out_buffer << "VERTICES " << tmp_points.size() << " " << tmp_points.size() * 2 << endl;
  for(m = 0; m < tmp_points.size(); m ++){
    out_buffer << "1 " << m << endl;
  }

  ofstream out_file(pm_filename.c_str());
  if (!out_file.is_open()){
    cerr << "ERROR: cannot open file " << pm_filename << endl;
    exit(1);
  }

  out_file << out_buffer.str();
  out_file.close();

}

// The first routine that the user needs to call to populate the surface data from file to 
// this data structure
void 
MySES_Raw :: init(string pm_dir_pdb, string pm_filename_pdb, string pm_filename_atom, string pm_filename_vert, string pm_filename_face){
  m_dir_pdb = pm_dir_pdb;
  m_filename_pdb = pm_filename_pdb;
  m_filename_xyzrn = pm_filename_atom;
  m_filename_vert = pm_filename_vert;
  m_filename_face = pm_filename_face;

  read_pdb_ca();

  //cout << "read_atom(...)..." << endl;
  read_xyzrn();

  //cout << "read_vert(...)..." << endl;
  read_vert();

  //cout << "read_face(...)..." << endl;
  read_face();

  read_patch();

  //cout << "get_bbox();" << endl;
  get_bbox();
}

// based on triangle info
void 
MySES_Raw :: read_patch() {
  size_t i;
  int patch_id, patch_type;

  m_patch_belt.clear();
  m_patch_concave.clear();
  m_patch_convex.clear();

  set<TID_VERT>   verts_empty;
  verts_empty.clear();

  for (i = 0; i < m_face_patch_id.size(); i ++){
    patch_type = m_face_type[i];
    patch_id = m_face_patch_id[i];
    if (patch_type == 1){
      if (m_patch_belt.find(patch_id) == m_patch_belt.end()){
        set<TID_VERT>   verts;
        verts.insert(m_face_triangle[i][0]);
        verts.insert(m_face_triangle[i][1]);
        verts.insert(m_face_triangle[i][2]);
        m_patch_belt[patch_id] = verts;
      }
      else{
        m_patch_belt[patch_id].insert(m_face_triangle[i][0]);
        m_patch_belt[patch_id].insert(m_face_triangle[i][1]);
        m_patch_belt[patch_id].insert(m_face_triangle[i][2]);
      }
    }
    else if (patch_type == 2){
      if (m_patch_concave.find(patch_id) == m_patch_concave.end()){
        set<TID_VERT>   verts;
        verts.insert(m_face_triangle[i][0]);
        verts.insert(m_face_triangle[i][1]);
        verts.insert(m_face_triangle[i][2]);
        m_patch_concave[patch_id] = verts;
      }
      else{
        m_patch_concave[patch_id].insert(m_face_triangle[i][0]);
        m_patch_concave[patch_id].insert(m_face_triangle[i][1]);
        m_patch_concave[patch_id].insert(m_face_triangle[i][2]);
      }

    }
    else if (patch_type == 3){
      if (m_patch_convex.find(patch_id) == m_patch_convex.end()){
        set<TID_VERT>   verts;
        verts.insert(m_face_triangle[i][0]);
        verts.insert(m_face_triangle[i][1]);
        verts.insert(m_face_triangle[i][2]);
        m_patch_convex[patch_id] = verts;
      }
      else{
        m_patch_convex[patch_id].insert(m_face_triangle[i][0]);
        m_patch_convex[patch_id].insert(m_face_triangle[i][1]);
        m_patch_convex[patch_id].insert(m_face_triangle[i][2]);
      }
    }
  }
} 


// Read face information from face file and fill out
// face related variables
void 
MySES_Raw :: read_face() {
        string in_filename = m_dir_pdb + m_filename_face; 
        ifstream in(in_filename.c_str());
        if (!in.is_open()){
            cerr << "ERROR: cannot open file " << in_filename << endl;
            exit(1);
        }

        char line[1024];
        line[1023] = '\0';

        TID_VERT id_v1, id_v2, id_v3;
        unsigned    patch_type, patch_id;
        size_t  i;
        size_t  num_atom;
        double  tmp_density;
        double  tmp_probe_radius;
        vector<TID_VERT>    tmp_triangle(3, 0);

        //discard the first 3 head lines
        in.getline(line, 1023);
        while (line[0] == '#' || strlen(line) == 0) in.getline(line, 1023);
        stringstream strbuf_0;
        strbuf_0 << line;
        strbuf_0 >> m_num_triangle >> num_atom >> tmp_density >> tmp_probe_radius;

        assert(num_atom == m_num_atom);
        assert(tmp_density == m_density);
        assert(tmp_probe_radius == m_probe_radius);

        m_face_triangle.clear();
        m_face_type.clear();
        m_face_patch_id.clear();

        m_face_triangle.reserve(m_num_triangle);
        m_face_type.reserve(m_num_triangle);
        m_face_patch_id.reserve(m_num_triangle);

        for (i = 0; i < m_num_triangle; i ++){       
            in.getline(line, 1023);
            stringstream strbuf;
            strbuf << line;
            strbuf >> id_v1 >> id_v2 >> id_v3 >> patch_type >> patch_id;

            tmp_triangle[0] = id_v1 - 1;
            tmp_triangle[1] = id_v2 - 1;
            tmp_triangle[2] = id_v3 - 1;

            //to fix the bug in .vert file, the info in .face is correct
            //m_vert_type[id_v1 - 1] = patch_type;
            //m_vert_type[id_v2 - 1] = patch_type;
            //m_vert_type[id_v3 - 1] = patch_type;

            m_face_triangle.push_back(tmp_triangle);
            m_face_type.push_back(4 - patch_type); //exchange 1 and 3
            m_face_patch_id.push_back(patch_id); 
        }
        in.close();
}

// Read vertex information from vert file and fill out related variables
void 
MySES_Raw :: read_vert() {
        string in_filename = m_dir_pdb + m_filename_vert; 
        ifstream in(in_filename.c_str());
        if (!in.is_open()){
            cerr << "ERROR: cannot open file " << in_filename << endl;
            exit(1);
        }
        char line[1024];
        line[1023] = '\0';
        double x, y, z, nx, ny, nz;
        int c7;
        TID_ATOM atomID;
        TID_VERT_TYPE vertType;
        string atomname;
        unsigned int    num_atom;
        vector<double>  xyz(3, 0);

        size_t i;

        in.getline(line, 1023);
        while (line[0] == '#') in.getline(line, 1023);
        stringstream strbuf0;
        strbuf0 << line;
        strbuf0 >> m_num_vertex >> num_atom >> m_density >> m_probe_radius; 

		if (num_atom != m_num_atom){
			cout << "num_atom = " << num_atom << endl;
			cout << "m_num_atom = " << m_num_atom << endl;
		}

		assert(num_atom == m_num_atom);

        m_vert_xyz.clear();         
        m_vert_normal.clear();  
        m_vert_c7.clear();      
        m_vert_atom_id.clear(); 
        m_vert_type.clear();    
        m_vert_atom_name.clear();

        m_vert_xyz.reserve(m_num_vertex);           
        m_vert_normal.reserve(m_num_vertex);  
        m_vert_c7.reserve(m_num_vertex);        
        m_vert_atom_id.reserve(m_num_vertex);   
        m_vert_type.reserve(m_num_vertex);  
        m_vert_atom_name.reserve(m_num_vertex);

        for (i = 0; i < m_num_vertex; i ++){
            in.getline(line, 1023);
            stringstream strbuf;
            strbuf << line;
            strbuf >> x >> y >> z 
                >> nx >> ny >> nz 
                >> c7 >> atomID >> vertType >> atomname;

            xyz[0] = x;
            xyz[1] = y;
            xyz[2] = z;
            m_vert_xyz.push_back(xyz);
            xyz[0] = nx;
            xyz[1] = ny;
            xyz[2] = nz;
            m_vert_normal.push_back(xyz);
            m_vert_c7.push_back(c7);
            m_vert_atom_id.push_back(atomID - 1);
            m_vert_type.push_back(vertType);
            m_vert_atom_name.push_back(atomname);
        }
        in.close();
}


// Read alpha carbon's co-ordinate from the pdb file
// and fill out alpah carbon related variables
void 
MySES_Raw :: read_pdb_ca() {
        string in_filename = m_dir_pdb + m_filename_pdb; 
        
        ifstream in(in_filename.c_str());
        if (!in.is_open()){
            cerr << "ERROR: cannot open file " << in_filename << ", unable to local ca from PDB file." << endl;
            exit(1);
        }
        
        //string out_filename = g_dir_output + m_filename_pdb + "_ca.txt";
        //ofstream out(out_filename.c_str());

        string line;
        vector<double>  tmp_xyz(3, 0);

        m_atom_xyz_ca.clear();
        m_atom_index_ca.clear();
        unsigned int  tmp_atom_index = 0;

        while(!in.eof()){
            getline(in, line);
            if (line.length() < 4) continue;
            if (line[0] == '#') continue;

            if (!((line[0] == 'A' || line[0] == 'a') 
                && (line[1] == 'T' || line[1] == 't') 
                && (line[2] == 'O' || line[2] == 'o') 
                && (line[3] == 'M' || line[3] == 'm')))
            {
                continue; 
            }
            tmp_atom_index ++;

            if (!((line[13] == 'C' || line[13] == 'c') && (line[14] == 'A' || line[14] == 'a'))) continue;

            tmp_xyz[0] = atof(line.substr(30, 8).c_str());   
            tmp_xyz[1] = atof(line.substr(38, 8).c_str());   
            tmp_xyz[2] = atof(line.substr(46, 8).c_str());   

            m_atom_xyz_ca.push_back(tmp_xyz);
            m_atom_index_ca.push_back(tmp_atom_index - 1);
            //out << line << endl;
        }    
        in.close();    
        //out.close();
}

// read xyzrn info from the xyzrn file and fill-out atom-related variables
void 
MySES_Raw :: read_xyzrn() {
        string in_filename = m_dir_pdb + m_filename_xyzrn; 
        ifstream in(in_filename.c_str());
        if (!in.is_open()){
            cerr << "ERROR: cannot open file " << in_filename << endl;
            exit(1);
        }

        char line[1024];
        line[1023] = '\0';
        double x, y, z, r;
        int c5;
        string atomname;  
        vector<double>  tmp_xyz(3, 0);

        m_atom_xyz.clear();
        m_atom_radius.clear();
        m_atom_c5.clear();
        m_atom_name.clear();

        while(!in.eof()){
            in.getline(line, 1023);
            if (strlen(line) == 0) continue;
            if (line[0] == '#') continue;

            stringstream strbuf;
            strbuf << line;
            strbuf >> x >> y >> z >> r >> c5 >> atomname;

            tmp_xyz[0] = x;
            tmp_xyz[1] = y;
            tmp_xyz[2] = z;

            m_atom_xyz.push_back(tmp_xyz);
            m_atom_radius.push_back(r);
            m_atom_c5.push_back(c5);
            m_atom_name.push_back(atomname);
        }    
        in.close();    

        m_num_atom = m_atom_xyz.size();
}

// find the bounding box and fill-out the variable m_min_x, m_max_y and so on
void 
MySES_Raw:: get_bbox() {
        double minx = 1000, miny = 1000, minz = 1000, maxx = -1000, maxy = -1000, maxz = -1000;
        size_t i;
        for (i = 0; i < m_vert_xyz.size(); i ++){
            if (minx > m_vert_xyz[i][0]) minx = m_vert_xyz[i][0];
            if (miny > m_vert_xyz[i][1]) miny = m_vert_xyz[i][1];
            if (minz > m_vert_xyz[i][2]) minz = m_vert_xyz[i][2];

            if (maxx < m_vert_xyz[i][0]) maxx = m_vert_xyz[i][0];
            if (maxy < m_vert_xyz[i][1]) maxy = m_vert_xyz[i][1];
            if (maxz < m_vert_xyz[i][2]) maxz = m_vert_xyz[i][2];
        }

        m_min_x = minx;
        m_max_x = maxx;
        m_min_y = miny;
        m_max_y = maxy;
        m_min_z = minz;
        m_max_z = maxz;

        //cout << "minx = " << minx << " miny = " << miny << " minz = " << minz << endl;
        //cout << "maxx = " << maxx << " maxy = " << maxy << " maxz = " << maxz << endl;
        //cout << "rangex = " << (maxx - minx) << " rangey = " << (maxy - miny) << " rangez = " << (maxz - minz) << endl;
}

// Given another MySES*, this routine defines the interface C-Alpha
void 
MySES :: populate_atom_interface_ca(MySES* pm_other_ses) {
        if (m_atom_interface_ca.size() > 0) return; // if current molecule already has C-Alpha filled, return

        vector<unsigned int>  tmp_local;

        size_t i;
        for (i = 0; i < m_atom_xyz_ca.size(); i ++){  // for each C-Alpha atom of current molecule
            pm_other_ses->query_local_atom(m_atom_xyz_ca[i], 10, tmp_local);  // find all the atoms of other_ses
	                                                                      // that are within 10A distance from
									      // this C-Alpha atom, store the result
									      // in tmp_local.
            if (tmp_local.size() > 0){
                m_atom_interface_ca.push_back(m_atom_xyz_ca[i]);              // since, there are atoms within 10A of 
		                                                              // this C-Alpha atom, it is interface C-Alpha
                m_atom_radius_ca.push_back(m_atom_radius[i]);
            }
        }
        cout << "ca# = " << m_atom_interface_ca.size() << endl;               // printing, How many of them?
}

// If you don't want to calculate interface C-Alpha all the time, you can load them from a file
// This routine load interface C-Alpha from a file
void 
MySES :: populate_atom_interface_ca() {
        m_atom_interface_ca.clear();
        m_atom_radius_ca.clear();

        string in_filename = g_dir_interface_ca + m_filename_pdb.substr(0, m_filename_pdb.size() - 4) + ".b.ca"; 
        ifstream in(in_filename.c_str());
        if (!in.is_open()){
            cerr << "ERROR: cannot open file " << in_filename << ", will calculate interface CA on fly." << endl;
            return;
        }

        string line;
        string  tmp_atom;
        int     tmp_atom_id;
        string  tmp_ca;
        string  tmp_atom_name;
        string  tmp_u;
        int     tmp_u_int;
        vector<double>  tmp_xyz(3, 0);

        m_atom_interface_ca.clear();
        m_atom_radius_ca.clear();

        while(!in.eof()){
            getline(in, line);
            if (line.length() < 4) continue;
            if (line[0] == '#') continue;

            stringstream strbuf;
            strbuf << line;
            strbuf >> tmp_atom 
                   >> tmp_atom_id 
                   >> tmp_ca 
                   >> tmp_atom_name
                   >> tmp_u
                   >> tmp_u_int
                   >> tmp_xyz[0]
                   >> tmp_xyz[1]
                   >> tmp_xyz[2];
            m_atom_interface_ca.push_back(tmp_xyz);
            m_atom_radius_ca.push_back(2);
        }    
        cout << "m_atom_interface_ca.size() = " << m_atom_interface_ca.size() << endl;
        in.close();    
}

// This routine populate the vertices from the variable m_vert_xyz, m_vert_type defined in MySES_Raw
void 
MySES :: populate_vert() {
        size_t i;
        MyVert* tmp_vert;

        m_vert.clear();
        m_vert.reserve(m_vert_xyz.size());

	// vertex id is given with i(1st arg of the constructor of MyVert)
	// type and xyz co-ordinates are coming from variables defined in MySES_Raw
        for (i = 0; i < m_vert_xyz.size(); i ++){
            tmp_vert = new MyVert(i, m_vert_type[i], m_vert_atom_id[i], m_vert_xyz[i][0], m_vert_xyz[i][1], m_vert_xyz[i][2]);
            m_vert.push_back(tmp_vert);
        }
}

// This routine populate the atoms from the variable m_atom_xyz[][], m_atom_radius, defined in MySES_Raw
void 
MySES :: populate_atom() {
        size_t i;
        MyAtom* tmp_atom;

        m_atom.clear();
        m_atom.reserve(m_atom_xyz.size());
        for (i = 0; i < m_atom_xyz.size(); i ++){
            tmp_atom = new MyAtom(i, m_atom_xyz[i][0], m_atom_xyz[i][1], m_atom_xyz[i][2], m_atom_radius[i]);
            m_atom.push_back(tmp_atom);
        }

	// Those atoms whise id, present in the m_atom_index_ca (define in MySES_Raw) are C-alpha atoms.
        for (i = 0; i < m_atom_index_ca.size(); i ++){
            m_atom[m_atom_index_ca[i]]->m_is_ca = true;
        }
}

// This routine populate the triangles from the variable m_face_triangle defined in MySES_Raw
void 
MySES :: populate_triangle() {   
        m_triangle = m_face_triangle;  
        vector<double>  v1(3, 0);  // a vertex point, <x,y, and z co-ordinates>
        vector<double>  v2(3, 0);  // a vertex point, <x,y, and z co-ordinates>
        vector<double>  v3(3, 0);  // a vertex point, <x,y, and z co-ordinates>
        double  tmp_area;

        TID_VERT    v1_id;  // a vertex id
        TID_VERT    v2_id;  // a vertex id
        TID_VERT    v3_id;  // a vertex id

        size_t i;

        for (i = 0; i < m_triangle.size(); i ++){
            v1_id = m_triangle[i][0];
            v2_id = m_triangle[i][1];
            v3_id = m_triangle[i][2];

            v1 = m_vert[v1_id]->point();  // read vertex co-ordinates are filling up
            v2 = m_vert[v2_id]->point();  // read vertex co-ordinates are filling up
            v3 = m_vert[v3_id]->point();  // read vertex co-ordinates are filling up

	    // Each vertex of the triangle gets 1/3 of the triangle area, a fair game
            tmp_area = compute_triangle_area(v1, v2, v3);
            m_vert[v1_id]->increase_area(tmp_area / 3);
            m_vert[v2_id]->increase_area(tmp_area / 3);
            m_vert[v3_id]->increase_area(tmp_area / 3);
        }
}

// A destructor
MySES :: ~MySES() {
        //m_index_vert->clear();
        delete m_index_vert;
        size_t i;
        for (i = 0; i < m_atom.size(); i ++){
            delete m_atom[i];
        }
        m_atom.clear();

        for (i = 0; i < m_vert.size(); i ++){
            delete m_vert[i];
        }
        m_vert.clear();
}

// I guess, this is indexing all the vertex point, using some fast data structure
void 
MySES :: build_index_vert() {
        vector<UINT_Key>        InputList;
        vector<double>          tmp_point(3, 0);
        size_t i;

        for (i = 0; i < m_vert.size(); i ++){
            tmp_point = m_vert[i]->point();
            InputList.push_back(UINT_Key(UINT_Pure_key(tmp_point[0], tmp_point[1], tmp_point[2]), i));
        }

        m_index_vert = new UINT_Range_tree_3_type(InputList.begin(), InputList.end());
}

// I guess, this is indexing all the center of the atoms, using some fast data structure
void 
MySES :: build_index_atom() {
        vector<UINT_Key>        InputList;
        vector<double>          tmp_point(3, 0);
        size_t i;

        for (i = 0; i < m_atom.size(); i ++){
            tmp_point = m_atom[i]->point();
            InputList.push_back(UINT_Key(UINT_Pure_key(tmp_point[0], tmp_point[1], tmp_point[2]), i));
        }

        m_index_atom = new UINT_Range_tree_3_type(InputList.begin(), InputList.end());
}

// For a given point(1st arg, pm_point) return all the vertex_id as a vector (3rd arg, pm_out),
// which are within the distance of pm_local_r, from pm_point
// CGAL, window_query is used. but window query works for rectangular region, so need to do a
// pruning at the end. (used for atom center)
void 
MySES :: query_local_atom(const vector<double>&  pm_point, double pm_local_r, vector<unsigned int>& pm_out) {
        pm_out.clear();
        pm_out.reserve(100);

        vector<UINT_Key> OutputList;
        double  tmp_local_dist = pm_local_r;

	// making windows
        UINT_Pure_key a = UINT_Pure_key(pm_point[0] - tmp_local_dist, pm_point[1] - tmp_local_dist, pm_point[2] - tmp_local_dist);
        UINT_Pure_key b = UINT_Pure_key(pm_point[0] + tmp_local_dist, pm_point[1] + tmp_local_dist, pm_point[2] + tmp_local_dist);
        UINT_Interval win = UINT_Interval(UINT_Key(a, 0), UINT_Key(b, 0));
	// runing queries
        m_index_atom->window_query(win, back_inserter(OutputList));

        if (OutputList.size() == 0) return;

        double tmp_dist;    
        size_t i;
	// doing pruning now.
        for (i = 0; i < OutputList.size(); i++){
            tmp_dist  = (pm_point[0] - OutputList[i].first.x()) * (pm_point[0] - OutputList[i].first.x());
            tmp_dist += (pm_point[1] - OutputList[i].first.y()) * (pm_point[1] - OutputList[i].first.y());
            tmp_dist += (pm_point[2] - OutputList[i].first.z()) * (pm_point[2] - OutputList[i].first.z());
            if (tmp_dist <= pm_local_r * pm_local_r){  // these points are within the local distance, hence copy in output
                pm_out.push_back(OutputList[i].second);
            }
        }

        return;
}

// For a given point(1st arg, pm_point) return all the vertex_id as a vector (3rd arg, pm_out),
// which are within the distance of pm_local_r, from pm_point
// CGAL, window_query is used. but window query works for rectangular region, so need to do a
// pruning at the end. (used for vertices)
void 
MySES :: query_local(const vector<double>&  pm_point, double pm_local_r, vector<unsigned int>& pm_out) {
        pm_out.clear();
        pm_out.reserve(3000);

        vector<UINT_Key> OutputList;
        double  tmp_local_dist = pm_local_r;

        UINT_Pure_key a = UINT_Pure_key(pm_point[0] - tmp_local_dist, pm_point[1] - tmp_local_dist, pm_point[2] - tmp_local_dist);
        UINT_Pure_key b = UINT_Pure_key(pm_point[0] + tmp_local_dist, pm_point[1] + tmp_local_dist, pm_point[2] + tmp_local_dist);
        UINT_Interval win = UINT_Interval(UINT_Key(a, 0), UINT_Key(b, 0));
        m_index_vert->window_query(win, back_inserter(OutputList));

        if (OutputList.size() == 0) return;

        double tmp_dist;    
        size_t i;
        for (i = 0; i < OutputList.size(); i++){
            tmp_dist  = (pm_point[0] - OutputList[i].first.x()) * (pm_point[0] - OutputList[i].first.x());
            tmp_dist += (pm_point[1] - OutputList[i].first.y()) * (pm_point[1] - OutputList[i].first.y());
            tmp_dist += (pm_point[2] - OutputList[i].first.z()) * (pm_point[2] - OutputList[i].first.z());
            if (tmp_dist <= pm_local_r * pm_local_r){
                pm_out.push_back(OutputList[i].second);
            }
        }

        return;
}
// Retrun the nearest neighbor of a point pm_point (used for vertices)
unsigned int 
MySES :: query_nn(const vector<double>& pm_point) {
        vector<UINT_Key> OutputList;

        double  tmp_local_dist = 1.0;
        UINT_Pure_key a = UINT_Pure_key(pm_point[0] - tmp_local_dist, pm_point[1] - tmp_local_dist, pm_point[2] - tmp_local_dist);
        UINT_Pure_key b = UINT_Pure_key(pm_point[0] + tmp_local_dist, pm_point[1] + tmp_local_dist, pm_point[2] + tmp_local_dist);
        UINT_Interval win = UINT_Interval(UINT_Key(a, 0), UINT_Key(b, 0));
        m_index_vert->window_query(win, back_inserter(OutputList));

	// incrementally increasing the local distances until we find at least one neighbor
        while (OutputList.size() == 0){
            tmp_local_dist += 0.2;
            a = UINT_Pure_key(pm_point[0] - tmp_local_dist, pm_point[1] - tmp_local_dist, pm_point[2] - tmp_local_dist);
            b = UINT_Pure_key(pm_point[0] + tmp_local_dist, pm_point[1] + tmp_local_dist, pm_point[2] + tmp_local_dist);
            win = UINT_Interval(UINT_Key(a, 0), UINT_Key(b, 0));
            m_index_vert->window_query(win, back_inserter(OutputList));
            //cerr << "ERROR: OutputList.size() == 0 in MyBall::query_nn(...)! " << endl;
            //cerr << "\twhen query (" << pm_point[0] << ", " << pm_point[1] << ", " << pm_point[2] << ")" << endl;
            //cerr << "\tm_nn_r = " << m_nn_r << endl;
            //double x = 10;
            //x = x / OutputList.size();
            ////exit(1);
        }

        vector<pair<double, unsigned int> > tmp_buffer;
        double tmp_dist;    
        size_t i;
        for (i = 0; i < OutputList.size(); i++){
            tmp_dist  = (pm_point[0] - OutputList[i].first.x()) * (pm_point[0] - OutputList[i].first.x());
            tmp_dist += (pm_point[1] - OutputList[i].first.y()) * (pm_point[1] - OutputList[i].first.y());
            tmp_dist += (pm_point[2] - OutputList[i].first.z()) * (pm_point[2] - OutputList[i].first.z());
            tmp_dist = sqrt(tmp_dist);
            tmp_buffer.push_back(pair<double, unsigned int>(tmp_dist, OutputList[i].second));
        }

        sort(tmp_buffer.begin(), tmp_buffer.end(), less_than_univ<double, unsigned int>());
        // Hasan: you just need the minimum, why r u sorting, just do find
        return tmp_buffer[0].second;
}



// Initializing routine for MySES. It first calls MySES_Raw::init, then it calls all the populate
// routines. Take all the input file names as argument
void 
MySES :: init(string pm_dir_pdb, string pm_filename_pdb, string pm_filename_atom, string pm_filename_vert, string pm_filename_face) {
        MySES_Raw::init(pm_dir_pdb, pm_filename_pdb, pm_filename_atom, pm_filename_vert, pm_filename_face);

        populate_atom();
        populate_vert();
        populate_triangle();
        build_index_vert();
        build_index_atom();
        //populate_atom_interface_ca();
}
