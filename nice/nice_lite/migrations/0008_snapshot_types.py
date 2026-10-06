# Stream particle sets name their snapshot kind: the 2D snapshots, typed 'snapshot' until
# now, become 'snapshot2D' beside multistate 3D's 'snapshot3D'. Rewrites the stored type in
# every job's particle_sets_stats.

from django.db import migrations


def _retype(apps, old, new):
    JobModel = apps.get_model('nice_lite', 'JobModel')
    for job in JobModel.objects.all().iterator():
        stats = job.particle_sets_stats
        particle_sets = stats.get('particle_sets') if isinstance(stats, dict) else None
        if not isinstance(particle_sets, list):
            continue
        changed = False
        for particle_set in particle_sets:
            if isinstance(particle_set, dict) and particle_set.get('type') == old:
                particle_set['type'] = new
                changed = True
        if changed:
            job.particle_sets_stats = stats
            job.save(update_fields=['particle_sets_stats'])


def snapshot_to_snapshot2D(apps, schema_editor):
    _retype(apps, 'snapshot', 'snapshot2D')


def snapshot2D_to_snapshot(apps, schema_editor):
    _retype(apps, 'snapshot2D', 'snapshot')


class Migration(migrations.Migration):

    dependencies = [
        ('nice_lite', '0007_rename_validate_projfile'),
    ]

    operations = [
        migrations.RunPython(snapshot_to_snapshot2D, snapshot2D_to_snapshot),
    ]
