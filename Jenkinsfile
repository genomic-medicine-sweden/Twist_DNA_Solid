def isPullRequest = false

pipeline {
  agent { label 'marvin-clone' }
  stages {
    stage('First stage') {
      steps {
        script {
          isPullRequest = env.CHANGE_ID != null
        }
        sh 'printenv'
      }
    }
    stage('HD832') {
      when {
          anyOf {
                 expression { isPullRequest == false && env.BRANCH_NAME == 'develop' }
          }
      }
      steps {
        sshagent(['jenkins']) {
          sh '''#!/bin/bash
                virtualenv venv -p python3.8
                source venv/bin/activate
                pip install -r requirements.txt
                mkdir -p HD832/slurm_out
                cp -r config HD832/
                cp .tests/release/units.tsv HD832/units.tsv
                cp .tests/release/samples_HD832.tsv HD832/samples.tsv
                cp .tests/release/resources.yaml HD832/config/resources.yaml
                cp .tests/release/config.yaml HD832/config/config.yaml

                module load singularity
                module load slurm-drmaa
                snakemake -s workflow/Snakefile --profile .tests/release/profile -d HD832 -w 60
             '''
         }
       }
    }
    stage('ALL') {
      when {
          anyOf {
                 expression { isPullRequest == false && env.BRANCH_NAME == 'master' }
          }
      }
      steps {
        sshagent(['jenkins']) {
          sh '''#!/bin/bash
                virtualenv venv -p python3.8
                source venv/bin/activate
                pip install -r requirements.txt
                mkdir -p ALL/slurm_out
                cp -r config ALL/
                cp .tests/release/units.tsv ALL/units.tsv
                cp .tests/release/samples_all.tsv ALL/samples.tsv
                cp .tests/release/resources.yaml ALL/config/resources.yaml
                cp .tests/release/config.yaml ALL/config/config.yaml

                module load singularity
                module load slurm-drmaa
                snakemake -s workflow/Snakefile --profile .tests/release/profile -d ALL -w 60
             '''
         }
       }
    }
  }
  post {
    always {
      cleanWs()

      dir("${env.WORKSPACE}@tmp") {
        deleteDir()
      }
      dir("${env.WORKSPACE}@script") {
        deleteDir()
      }
      dir("${env.WORKSPACE}@script@tmp") {
        deleteDir()
      }
    }
  }
}
