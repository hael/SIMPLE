if(${NICE})

    find_package(Python3 3.10)

    add_custom_target(nice_src ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different ${CMAKE_CURRENT_SOURCE_DIR}/example_templates ${CMAKE_CURRENT_BINARY_DIR}/example_templates
        COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different ${CMAKE_CURRENT_SOURCE_DIR}/example_systemd   ${CMAKE_CURRENT_BINARY_DIR}/example_systemd
        COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different ${CMAKE_CURRENT_SOURCE_DIR}/example_proxy     ${CMAKE_CURRENT_BINARY_DIR}/example_proxy
        COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different ${CMAKE_CURRENT_SOURCE_DIR}/nice_project      ${CMAKE_CURRENT_BINARY_DIR}/nice_project
        COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different ${CMAKE_CURRENT_SOURCE_DIR}/nice_lite         ${CMAKE_CURRENT_BINARY_DIR}/nice_lite
        COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different ${CMAKE_CURRENT_SOURCE_DIR}/test_display      ${CMAKE_CURRENT_BINARY_DIR}/test_display
        COMMAND ${CMAKE_COMMAND} -E copy                        ${CMAKE_CURRENT_SOURCE_DIR}/manage.py         ${CMAKE_CURRENT_BINARY_DIR}/manage.py
        COMMAND ${CMAKE_COMMAND} -D NICE_FQDN=${NICE_FQDN} -D SIMPLE_PATH=${CMAKE_BINARY_DIR} -D SRC_PATH=${CMAKE_CURRENT_SOURCE_DIR} -P "${CMAKE_CURRENT_LIST_DIR}/substitute_nice_src.cmake"
    )

    if(${APPLE})
        add_custom_target(nice_venv ALL
            COMMAND ${Python3_EXECUTABLE} -m venv ${CMAKE_CURRENT_BINARY_DIR}/venv
            COMMAND ./venv/bin/python3 -m pip install --upgrade pip
            COMMAND ./venv/bin/pip3 install --isolated --require-virtualenv --no-cache-dir --no-input django==5.2.5 pytailwindcss==0.2.0 dirsync==2.2.6 
        )
    else()
        add_custom_target(nice_venv ALL
            COMMAND ${Python3_EXECUTABLE} -m venv ${CMAKE_CURRENT_BINARY_DIR}/venv
            COMMAND ./venv/bin/python3 -m pip install --upgrade pip
            COMMAND ./venv/bin/pip3 install --isolated --require-virtualenv --no-cache-dir --no-input django==5.2.5 pytailwindcss==0.2.0 dirsync==2.2.6 pysqlite3-binary==0.5.4 gunicorn==23.0.0
        )
    endif()

    add_custom_target(nice_build ALL
        DEPENDS nice_src
        DEPENDS nice_venv
        COMMAND ./venv/bin/tailwindcss -i nice_lite/templates/tailwind.css -o nice_lite/static/nice_lite/nice_lite.css --minify
        COMMAND ./venv/bin/python3 ${CMAKE_CURRENT_BINARY_DIR}/manage.py makemigrations
        COMMAND ./venv/bin/python3 ${CMAKE_CURRENT_BINARY_DIR}/manage.py migrate
        COMMAND ./venv/bin/python3 ${CMAKE_CURRENT_BINARY_DIR}/manage.py collectstatic --no-input
    )

endif()
