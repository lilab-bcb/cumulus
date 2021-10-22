Import workflows
------------------------------------

We list here 2 ways to import and execute Cumulus workflows

Import from Dockstore to Terra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Cumulus workflows are hosted on Dockstore_ under organization of Broad Institute of MIT and Harvard. Select "Cumulus Workflows" in Dockstore. 

.. image:: images/select_workflow.png
   :scale: 60 %
   :align: center

2. For purpose of illustration, we will describe process of importing Demultiplexing workflow in Terra. Click on View button beside "github.com/klarman-cell-observatory/cumulus/Demultiplexing".

.. image:: images/demult_workflow.png
   :scale: 60 %
   :align: center

3. Switch version using the Versions tab. In "Git Reference column" select the appropriate version to import into Terra. 

.. image:: images/version_selection.png
   :scale: 60 %
   :align: center

4. Then Launch with Terra (illustrated below).

.. image:: images/launch_terra.png
   :scale: 60 %
   :align: center

5. Type in an appropriate "Workflow Name" & Select a workspace on Terra where you have access to execute workflows.    

Command line usage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another method of launching workflow is by using command line tool Altocumulus_ . After installing Altocumulus use command line below ::

    alto cromwell run -s 10.10.10.10 -p 3000 -m broadinstitute:cumulus:Demultiplexing:master \
                      -i cumulus_inputs.json --no-ssl-verify

where:

      - -m specifies which WDL workflow to use. It uses the Dockstore-style specification. In the example above, Cumulus workflow will be used, and in specific it's the master branch version to be used. If the version is left blank, altocumulus will use its default version on Dockstore.
      - -s specifies the internal IP address of server.
      - -p specifies the port.
      - -i specifies the local workflow input JSON file.
      - --no-ssl-verify: This option is needed when in Roche intranet.


.. _Dockstore: https://dockstore.org/organizations/BroadInstitute
.. _Altocumulus: https://github.com/klarman-cell-observatory/altocumulus/tree/master/alto/commands
