from PyMEMENTO import MEMENTO

chain_ids = []
for letter in ["A","B","C","D","E"]:
    chain_ids += [letter for i in range(390)]
residue_ids = (list(range(1,321))+list(range(409,479)))*5

model = MEMENTO(
    "testrun/",
    "testing_alpha7/desensitized/protein_densensitized.gro",
    "testing_alpha7/open/protein_lipids_open.gro",
    multiple_chains=chain_ids,
    last_step_performed='processing',
    lipid="resname PA or resname PC or resname OL"
)
# Perform morphing and modelling
#model.morph(24)
#model.make_models(50)
#model.find_best_path()
#model.process_models(residue_ids, caps=False)
model.prepare_boxes(template_folder="testing_alpha7/template_alpha7", mdrun_flags={'ntomp':6,'ntmpi':1,'gpu_id':0},embedding_starting_scale=1.3, embedding_end_scale=1, embedding_steps=5)
model.solvate_boxes(ion_concentration=0.15)
model.minimize_boxes(mdrun_flags={'ntomp':6,'ntmpi':1,'gpu_id':0, 'pin':'on', 'pinoffset':0}, grompp_flags={'maxwarn':1})