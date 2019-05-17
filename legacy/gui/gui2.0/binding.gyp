{
  "targets": [
    {
      "include_dirs": [
        "<!(node -e \"require('nan')\")"
      ],
      "target_name": "MRChandler",
      "sources": [ "src/addons/MRChandler.cpp" ]
    }
  ]
}
