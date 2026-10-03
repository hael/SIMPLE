# 2026-10-03 program renames: abinitio2D/abinitio3D became solve2D/solve3D (abinitio2D_stream
# became pool2D) and cluster2D became refine2D.
# Renames the multistate stream-stage fields and rewrites stored program names.

from django.db import migrations

# legacy program name -> current program name (mirrors simple_ui_legacy_names.f90)
LEGACY_PROGS = {
    'abinitio2D':                 'solve2D',
    'abinitio2D_chunks':          'solve2D_chunks',
    'abinitio2D_stream':          'pool2D',
    'abinitio3D':                 'solve3D',
    'abinitio3D_cavgs':           'solve3D_cavgs',
    'abinitio3D_addon':           'solve3D_addon',
    'abinitio3D_addon_snapshots': 'solve3D_addon_snapshots',
    'abinitio3D_nano':            'solve3D_nano',
    'abinitio3D_stream':          'solve3D_stream',
    'cluster2D':                  'refine2D',
    'cluster2D_distr':            'refine2D_distr',
    'cluster2D_nano':             'refine2D_nano',
}


def rename_progs(apps, schema_editor):
    JobModel = apps.get_model('nice_lite', 'JobModel')
    for old, new in LEGACY_PROGS.items():
        JobModel.objects.filter(prog=old).update(prog=new)


def restore_progs(apps, schema_editor):
    JobModel = apps.get_model('nice_lite', 'JobModel')
    for old, new in LEGACY_PROGS.items():
        JobModel.objects.filter(prog=new).update(prog=old)


class Migration(migrations.Migration):

    dependencies = [
        ('nice_lite', '0005_jobmodel_abinitio3d_multistate_heartbeat_and_more'),
    ]

    operations = [
        migrations.RenameField(
            model_name='jobmodel',
            old_name='abinitio3D_multistate_heartbeat',
            new_name='solve3D_multistate_heartbeat',
        ),
        migrations.RenameField(
            model_name='jobmodel',
            old_name='abinitio3D_multistate_stats',
            new_name='solve3D_multistate_stats',
        ),
        migrations.RenameField(
            model_name='jobmodel',
            old_name='abinitio3D_multistate_status',
            new_name='solve3D_multistate_status',
        ),
        migrations.RenameField(
            model_name='jobmodel',
            old_name='abinitio3D_multistate_update',
            new_name='solve3D_multistate_update',
        ),
        migrations.RunPython(rename_progs, restore_progs),
    ]
