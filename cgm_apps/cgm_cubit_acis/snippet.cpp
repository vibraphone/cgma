#include "RefEntityName.hpp"
#include "RefGroup.hpp"
#include "RefFace.hpp"
#include "CubitString.hpp"
#include <iostream>

int snippet()
{

  RefEntity* this_grp = RefEntityName::instance()->get_refentity(CubitString("spec_reflect"));
  RefGroup *this_group = dynamic_cast<RefGroup*>(this_grp);
  if (NULL == this_group) return 1;
  DLIList<RefEntity*> grp_ents;
  this_group->expand_group(grp_ents);
  for (int i = grp_ents.size(); i > 0; i--) {
    RefEntity *this_ent = grp_ents.get_and_step();
      // don't really need this_face unless you want to check its type
    RefFace *this_face = dynamic_cast<RefFace*>(this_ent);
    if (NULL == this_face) continue;
    int this_id = this_ent->id();
    std::cout << "Surface " << this_id << " found." << std::endl;
  }

  return 0;
}
