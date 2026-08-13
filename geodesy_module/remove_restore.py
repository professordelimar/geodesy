"""Fluxo remove-restore para modelagem regional."""
import numpy as np
def remove_global_model(observed,global_model): return np.asarray(observed,dtype=float)-np.asarray(global_model,dtype=float)
def restore_global_model(residual,global_model): return np.asarray(residual,dtype=float)+np.asarray(global_model,dtype=float)
def remove_terrain_effect(observed,terrain_effect): return np.asarray(observed,dtype=float)-np.asarray(terrain_effect,dtype=float)
def restore_terrain_effect(residual,terrain_effect): return np.asarray(residual,dtype=float)+np.asarray(terrain_effect,dtype=float)
def remove_restore_workflow(observed,global_model=None,terrain_effect=None,modeling_function=None,**kwargs):
    residual=np.asarray(observed,dtype=float).copy()
    if global_model is not None: residual=remove_global_model(residual,global_model)
    if terrain_effect is not None: residual=remove_terrain_effect(residual,terrain_effect)
    modeled_residual=modeling_function(residual,**kwargs) if modeling_function is not None else residual
    final=modeled_residual
    if terrain_effect is not None: final=restore_terrain_effect(final,terrain_effect)
    if global_model is not None: final=restore_global_model(final,global_model)
    return final,residual,modeled_residual
