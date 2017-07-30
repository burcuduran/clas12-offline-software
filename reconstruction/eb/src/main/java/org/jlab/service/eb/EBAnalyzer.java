package org.jlab.service.eb;

import static java.lang.Math.abs;
import static java.lang.Math.pow;
import java.util.HashMap;
import org.jlab.clas.detector.DetectorEvent;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.detector.base.DetectorType;


/**
 *
 * @author gavalian
 * @author jnewton
 */
public class EBAnalyzer {
    
    private int[] pidPositives = new int[]{-11,  211, 321, 2212};
    private int[] pidNegatives = new int[]{ 11, -211,-321,-2212};
    
    public EBAnalyzer(){}
    
    public void processEvent(DetectorEvent event) {

        if (event.getParticles().size()==0) return;

        DetectorParticle trigger = event.getParticle(0);

        // assign start time based on trigger particle:
        if (trigger.getPid()==11 || trigger.getPid()==-11) {

            trigger.setBeta(1.0);
            trigger.setMass(0.0005);

            double time = 0.0;
            double path = 0.0;

            // FTOF Panel 1B:
            if(trigger.hasHit(DetectorType.FTOF, 2)){
                time = trigger.getTime(DetectorType.FTOF, 2);
                path = trigger.getPathLength(DetectorType.FTOF, 2);
            }
            // FTOF Panel 1A:
            else if(trigger.hasHit(DetectorType.FTOF, 1)){
                time = trigger.getTime(DetectorType.FTOF, 1);
                path = trigger.getPathLength(DetectorType.FTOF, 1);
            }

            double tof = path/EBConstants.SPEED_OF_LIGHT;

            double start_time = time - tof;

            double deltatr = - start_time + event.getEventHeader().getRfTime()
                    + (EBConstants.RF_LARGE_INTEGER+0.5)*EBConstants.RF_BUCKET_LENGTH + EBConstants.RF_OFFSET;

            double rfcorr = deltatr%EBConstants.RF_BUCKET_LENGTH - EBConstants.RF_BUCKET_LENGTH/2;

            event.getEventHeader().setStartTime(start_time + rfcorr);

            //double deltatr = - start_time + event.getEventHeader().getRfTime() /* - (trigger.vertex().z()
            //        - (EBConstants.TARGET_POSITION))/(EBConstants.SPEED_OF_LIGHT)*/
            //        + (EBConstants.RF_LARGE_INTEGER+0.5)*EBConstants.RF_BUCKET_LENGTH + EBConstants.RF_OFFSET;

            //System.out.println(event.getEventHeader().getRfTime() - start_time);
            //System.out.println(rfcorr + " " + (124.25- time + tof));
            //System.out.println(" TIME = " + tof + "  time from TOF = " + time);
            //System.out.println(" PATH = " + path + " " );
            //System.out.println(" SET START TIME = " + start_time + "  ACTUAL TIME = " + event.getEventHeader().getStartTime());
            //System.out.println(start_time - event.getEventHeader().getRfTime());
        }
        else if (trigger.getPid()==0 || trigger.getPid()==22) {
            event.getEventHeader().setStartTime(124.25);
        }

        // identify non-trigger particlef:
        this.assignMasses(event);
        this.assignPids(event);
    }
    
    
    public void assignMasses(DetectorEvent event){
        
        // start at second particle, because
        // first is trigger particle with PID assigned elsewhere:
        for(int i = 1; i < event.getParticles().size(); i++) {
            DetectorParticle p = event.getParticle(i);
            double start_time  = event.getEventHeader().getStartTime();
            double beta = 0.0;
            double mass = 0.0;

            // FTOF Panel 1B:
            if(p.hasHit(DetectorType.FTOF, 2)){
                beta = p.getBeta(DetectorType.FTOF, 2,start_time);
                //mass = p.getMass2(DetectorType.FTOF, 2,start_time);
                p.setBeta(beta);
            }

            // FTOF Panel 1A:
            else if(p.hasHit(DetectorType.FTOF, 1)){
                beta = p.getBeta(DetectorType.FTOF,1, start_time);
                //mass = p.getMass2(DetectorType.FTOF,1, start_time);
                p.setBeta(beta);
            }

            // CTOF:
            else if(p.hasHit(DetectorType.CTOF, 0)){
                beta = p.getBeta(DetectorType.CTOF ,start_time);
                //mass = p.getMass2(DetectorType.CTOF,start_time);
                p.setBeta(beta);
            }
        }
    }
    
    public void assignPids(DetectorEvent event) {

        PIDHypothesis pidHyp = new PIDHypothesis();
        pidHyp.setEvent(event);

        // start at second particle, because first particle is the
        // trigger particle with PID assigned elsewhere:
        for (int i = 1; i < event.getParticles().size(); i++) {

            final int pid=pidHyp.getSimplePid(event.getParticle(i));
            pidHyp.finalizePID(event.getParticle(i),pid);
        }

        /*
        for (int i = 1; i < event.getParticles().size(); i++){

            DetectorParticle p = event.getParticle(i);

            // Neutrals:
            if(p.getCharge()==0) {
                break;
            }

            // Positives:
            else if(p.getCharge()>0){
                for (int aPid : pidPositives) {
                    pidHyp.setEvent(event);
                    pidHyp.PIDMatch(p, aPid);
                    //pidHyp.PIDQuality(p,aPid,event);
                }
            }

            // Negatives:
            else {
                for (int aPid : pidNegatives) {
                    pidHyp.setEvent(event);
                    pidHyp.PIDMatch(p, aPid);
                    //pidHyp.PIDQuality(p, aPid,event);
                }
            }
        }
        */
    }


    public class PIDHypothesis {

        private int theoryPID = -1;
        private double PIDquality = 0.0;
        private DetectorEvent event;

        public PIDHypothesis() {
        }

        public void setEvent(DetectorEvent e) {
            event = e;
        }

        public void PIDMatch(DetectorParticle p, int pid) {

            //double vertex_index = optimalVertexTime(p);
            //boolean vertexCheck = (abs(pid)==211 && vertex_index==1 && p.getBeta()>0.0) ||
            //        (abs(pid)==2212 && vertex_index==0 && p.getBeta()>0.0) ||
            //        (abs(pid)==321 && vertex_index==2 && p.getBeta()>0.0);

            boolean vertexCheck = (abs(pid) == bestPidBasedOnTiming(p) && p.getBeta() > 0.0);

            boolean sfCheck = p.getEnergyFraction(DetectorType.EC) > EBConstants.ECAL_SAMPLINGFRACTION_CUT;
            boolean htccSignalCheck = p.getNphe(DetectorType.HTCC) > EBConstants.HTCC_NPHE_CUT;
            boolean ltccSignalCheck = p.getNphe(DetectorType.LTCC) > EBConstants.LTCC_NPHE_CUT;
            boolean htccPionThreshold = p.vector().mag() > EBConstants.HTCC_PION_THRESHOLD;
            boolean ltccPionThreshold = p.vector().mag() < EBConstants.LTCC_UPPER_PION_THRESHOLD
                    && p.vector().mag() > EBConstants.LTCC_LOWER_PION_THRESHOLD;

            switch (abs(pid)) {
                case 11:
                    if (htccSignalCheck && sfCheck) {
                        this.finalizePID(p, pid);
                        break;
                    }
                case 211:
                    if (htccSignalCheck && !sfCheck && htccPionThreshold) {
                        this.finalizePID(p, pid);
                        break;
                    }
                    if (ltccSignalCheck && !sfCheck && ltccPionThreshold) {
                        this.finalizePID(p, pid);
                        break;
                    }
                case 321:
                    if (vertexCheck && !sfCheck && !htccSignalCheck) {
                        this.finalizePID(p, pid);
                        break;
                    }
                case 2212:
                    if (vertexCheck && !sfCheck && !htccSignalCheck) {
                        this.finalizePID(p, pid);
                        break;
                    }
            }
        }

        public int getSimplePid(DetectorParticle p) {

            boolean sfCheck = p.getEnergyFraction(DetectorType.EC) > EBConstants.ECAL_SAMPLINGFRACTION_CUT;
            boolean htccSignalCheck = p.getNphe(DetectorType.HTCC) > EBConstants.HTCC_NPHE_CUT;
            boolean ltccSignalCheck = p.getNphe(DetectorType.LTCC) > EBConstants.LTCC_NPHE_CUT;
            boolean htccPionThreshold = p.vector().mag() > EBConstants.HTCC_PION_THRESHOLD;
            //boolean ltccPionThreshold = p.vector().mag() > EBConstants.LTCC_LOWER_PION_THRESHOLD;
            boolean ltccKaonThreshold = p.vector().mag() > EBConstants.LTCC_UPPER_PION_THRESHOLD;

            final int qq = p.getCharge();

            int simplePid = 0;

            // not ready for neutrals:
            if (qq == 0) simplePid = 0;

            // charged:
            else {

                // electron/positron requirements:
                if (htccSignalCheck && sfCheck) simplePid = qq * -11;

                // else it's a hadron:
                else {
                    final int timingPid = bestPidBasedOnTiming(p);

                    switch (timingPid) {
                        case 211:
                            simplePid = qq * 211;
                            break;
                        case 321:
                            if (ltccSignalCheck && !ltccKaonThreshold) simplePid = qq * 211;
                            else if (htccSignalCheck && htccPionThreshold) simplePid = qq * 211;
                            else simplePid = qq * 321;
                            break;
                        case 2212:
                            if (ltccSignalCheck && ltccKaonThreshold) simplePid = qq * 321;
                            else if (p.getCharge() > 0) simplePid = 2212;
                            else simplePid = qq * 321;
                            break;
                        default:
                            throw new RuntimeException("EBAnalyzer:bestPid:  Unknown pid:  " + timingPid);
                    }
                }
            }
            return simplePid;
        }

/*
        public int optimalVertexTime(DetectorParticle p) {
            int vertex_index = 0;
            HashMap<Integer,Double> vertexDiffs = new HashMap<Integer,Double>();
            double event_start_time = event.getEventHeader().getStartTime();

            // FTOF Panel 1B:
            if(p.hasHit(DetectorType.FTOF,2)) {
                vertexDiffs.put(0,abs(p.getVertexTime(DetectorType.FTOF, 2, 2212)-event_start_time));
                vertexDiffs.put(1,abs(p.getVertexTime(DetectorType.FTOF, 2, 211)-event_start_time));
                vertexDiffs.put(2,abs(p.getVertexTime(DetectorType.FTOF, 2, 321)-event_start_time));
            }

            // FTOF Panel 1A:
            else if(p.hasHit(DetectorType.FTOF,1)) {
                vertexDiffs.put(0,abs(p.getVertexTime(DetectorType.FTOF, 1, 2212)-event_start_time));
                vertexDiffs.put(1,abs(p.getVertexTime(DetectorType.FTOF, 1, 211)-event_start_time));
                vertexDiffs.put(2,abs(p.getVertexTime(DetectorType.FTOF, 1, 321)-event_start_time));
            }

            if(vertexDiffs.size()>0) {
                double min = vertexDiffs.get(0);
                for (int i = 1; i <= 2; i++) {
                    if (vertexDiffs.get(i) < min) {
                        min = vertexDiffs.get(i);
                        vertex_index = i;
                    }
                }
            }
            return vertex_index;
        }
*/

        public int bestPidBasedOnTiming(DetectorParticle p) {

            int bestPid=0;

            // Prefer FTOF Panel 1B (2) over 1A (1):
            int iPanel=0;
            if     (p.hasHit(DetectorType.FTOF,2)) iPanel=2;
            else if(p.hasHit(DetectorType.FTOF,1)) iPanel=1;

            if (iPanel>0) {
                double minTimeDiff=Double.MAX_VALUE;
                double startTime = event.getEventHeader().getStartTime();
                for (int aPid : new int [] {2212,211,321}) {
                    double timeDiff=abs(p.getVertexTime(DetectorType.FTOF, iPanel, aPid)-startTime);
                    if (timeDiff < minTimeDiff) {
                        bestPid = aPid;
                        minTimeDiff = timeDiff;
                    }
                }
            }

            return bestPid;
        }

        public double PIDQuality(DetectorParticle p, int pid, DetectorEvent event) {
            double delta_t = abs(p.getVertexTime(DetectorType.FTOF, 2, pid)-event.getEventHeader().getStartTime());
            double sigma = 0.08;
            return pow((delta_t/sigma),2);
        }

        public void finalizePID(DetectorParticle p, int pid) {
            p.setPid(pid);
            theoryPID = pid;
            PIDquality = this.PIDQuality(p, pid, event);
            p.setPidQuality(PIDquality);
        }
    }

}



