#pragma once
#include "Constants.h"
#include "Ensemble.h"
#include "Grid.h"
#include "AnalyticalField.h"
#include "Pusher.h"
#include "Vectors.h"

#include <stdint.h>
#include <omp.h>
#include <random>

using namespace constants;
namespace pfc
{
    template <class TGrid>  // may be AnalyticalField or any Grid type
    class Scalar_Fast_QED : public ParticlePusher
    {
    public:

        Scalar_Fast_QED()
        {
            g_emis = (double*)int_g_emis;
            g_pair = (double*)int_g_pair;

            SchwingerField = sqr(Constants<FP>::electronMass() * Constants<FP>::lightVelocity())
                * Constants<FP>::lightVelocity() / (-Constants<FP>::electronCharge() * Constants<FP>::planck());

            preFactor = sqr(Constants<FP>::electronCharge()) * Constants<FP>::electronMass()
                * Constants<FP>::lightVelocity() / sqr(Constants<FP>::planck());

            coeffPhoton_probability = 1.0;
            coeffPair_probability = 1.0;

            distribution = std::uniform_real_distribution<FP>(0.0, 1.0);
            int max_threads;

#ifdef __USE_OMP__
            max_threads = omp_get_max_threads();
#else
            max_threads = 1;
#endif
            AvalanchePhotons.resize(max_threads);
            AvalancheParticles.resize(max_threads);
            afterAvalanchePhotons.resize(max_threads);
            afterAvalancheParticles.resize(max_threads);
            timeAvalanchePhotons.resize(max_threads);
            timeAvalancheParticles.resize(max_threads);
        }

        void disable_photon_emission() 
        {
            this->coeffPhoton_probability = 0.0;
        }

        void enable_photon_emission()
        {
            this->coeffPhoton_probability = 1.0;
        }

        void disable_pair_production()
        {
            this->coeffPair_probability = 0.0;
        }

        void enable_pair_production()
        {
            this->coeffPair_probability = 1.0;
        }

        void processParticles(Ensemble3d* particles, TGrid* grid, FP timeStep)
        {
            int max_threads;
#ifdef __USE_OMP__
            max_threads = omp_get_max_threads();
#else
            max_threads = 1;
#endif

            for (int th = 0; th < max_threads; th++)
            {
                AvalanchePhotons[th].clear();
                AvalancheParticles[th].clear();
                afterAvalanchePhotons[th].clear();
                afterAvalancheParticles[th].clear();
            }

            if ((*particles)[Photon].size())
                HandlePhotons((*particles)[Photon], grid, timeStep);
            if ((*particles)[Electron].size())
                HandleParticles((*particles)[Electron], grid, timeStep);
            if ((*particles)[Positron].size())
                HandleParticles((*particles)[Positron], grid, timeStep);

            for (int th = 0; th < max_threads; th++)
            {
                for (int ind = 0; ind < afterAvalanchePhotons[th].size(); ind++)
                {
                    particles->addParticle(afterAvalanchePhotons[th][ind]);
                }
                for (int ind = 0; ind < afterAvalancheParticles[th].size(); ind++)
                {
                    particles->addParticle(afterAvalancheParticles[th][ind]);
                }
            }
        }

        void Boris(Particle3d&& particle, const FP3& e, const FP3& b, FP timeStep)
        {
            FP eCoeff = timeStep * particle.getCharge() / (2 * particle.getMass() * Constants<FP>::lightVelocity());
            FP3 eMomentum = e * eCoeff;
            FP3 um = particle.getP() + eMomentum;
            FP3 t = b * eCoeff / sqrt((FP)1 + um.norm2());
            FP3 uprime = um + cross(um, t);
            FP3 s = t * (FP)2 / ((FP)1 + t.norm2());
            particle.setP(eMomentum + um + cross(uprime, s));
            particle.setPosition(particle.getPosition() + timeStep * particle.getVelocity());
        }

        void Boris(ParticleProxy3d&& particle, const FP3& e, const FP3& b, FP timeStep)
        {
            FP eCoeff = timeStep * particle.getCharge() / (2 * particle.getMass() * Constants<FP>::lightVelocity());
            FP3 eMomentum = e * eCoeff;
            FP3 um = particle.getP() + eMomentum;
            FP3 t = b * eCoeff / sqrt((FP)1 + um.norm2());
            FP3 uprime = um + cross(um, t);
            FP3 s = t * (FP)2 / ((FP)1 + t.norm2());
            particle.setP(eMomentum + um + cross(uprime, s));
            particle.setPosition(particle.getPosition() + timeStep * particle.getVelocity());
        }

        void HandlePhotons(ParticleArray3d& particles, TGrid* grid, FP timeStep)
        {
#pragma omp parallel for
                for (int i = 0; i < particles.size(); i++)
                {
                    int thread_id;
#ifdef __USE_OMP__
                    thread_id = omp_get_thread_num();
#else
                    thread_id = 0;
#endif
                    FP3 pPos = particles[i].getPosition();
                    FP3 e, b;

                    e = grid->getE(pPos);
                    b = grid->getB(pPos);

                    AvalancheParticles[thread_id].clear();
                    AvalanchePhotons[thread_id].clear();
                    timeAvalancheParticles[thread_id].clear();
                    timeAvalanchePhotons[thread_id].clear();
                    AvalanchePhotons[thread_id].push_back(particles[i]);
                    timeAvalanchePhotons[thread_id].push_back((FP)0.0);

                    RunAvalanche(e, b, timeStep);

                    for (int k = 0; k != AvalanchePhotons[thread_id].size(); k++)
                        if (AvalanchePhotons[thread_id][k].getGamma() != (FP)1.0)
                            afterAvalanchePhotons[thread_id].push_back(AvalanchePhotons[thread_id][k]);


                    for (int k = 0; k != AvalancheParticles[thread_id].size(); k++)
                        afterAvalancheParticles[thread_id].push_back(AvalancheParticles[thread_id][k]);
                }

                particles.clear();
        }

        void HandleParticles(ParticleArray3d& particles, TGrid* grid, FP timeStep)
        {
            
#pragma omp parallel for
            for (int i = 0; i < particles.size(); i++)
            {
                int thread_id;
#ifdef __USE_OMP__
                thread_id = omp_get_thread_num();
#else
                thread_id = 0;
#endif
                FP3 pPos = particles[i].getPosition();
                FP3 e, b;

                e = grid->getE(pPos);
                b = grid->getB(pPos);
                
                AvalancheParticles[thread_id].clear();
                AvalanchePhotons[thread_id].clear();
                timeAvalancheParticles[thread_id].clear();
                timeAvalanchePhotons[thread_id].clear();
                AvalancheParticles[thread_id].push_back(particles[i]);
                timeAvalancheParticles[thread_id].push_back((FP)0.0);

                RunAvalanche(e, b, timeStep);

                for (int k = 0; k != AvalanchePhotons[thread_id].size(); k++)
                    if(AvalanchePhotons[thread_id][k].getGamma() != (FP)1.0)
                        afterAvalanchePhotons[thread_id].push_back(AvalanchePhotons[thread_id][k]);

                particles[i].setMomentum(AvalancheParticles[thread_id][0].getMomentum());
                particles[i].setPosition(AvalancheParticles[thread_id][0].getPosition());

                for (int k = 1; k != AvalancheParticles[thread_id].size(); k++)
                    afterAvalancheParticles[thread_id].push_back(AvalancheParticles[thread_id][k]);
            }
        }
        
        void RunAvalanche(const FP3& E, const FP3& B, double timeStep)
        {
            int thread_id;
#ifdef __USE_OMP__
            thread_id = omp_get_thread_num();
#else
            thread_id = 0;
#endif
            vector<Particle3d>& AvalancheParticles = this->AvalancheParticles[thread_id];
            vector<Particle3d>& AvalanchePhotons = this->AvalanchePhotons[thread_id];

            vector<FP>& timeAvalancheParticles = this->timeAvalancheParticles[thread_id];
            vector<FP>& timeAvalanchePhotons = this->timeAvalanchePhotons[thread_id];
            
            int countParticles = 0;
            int countPhotons = 0;

            while (countParticles != AvalancheParticles.size()
                || countPhotons != AvalanchePhotons.size())
            {
                for (int k = countParticles; k != AvalancheParticles.size(); k++)
                {
                    oneParticleStep(AvalancheParticles[k], E, B, timeAvalancheParticles[k], timeStep);
                    countParticles++;
                }
                for (int k = countPhotons; k != AvalanchePhotons.size(); k++)
                {
                    onePhotonStep(AvalanchePhotons[k], E, B, timeAvalanchePhotons[k], timeStep);
                    countPhotons++;
                }
            }
        }

        void oneParticleStep(Particle3d& particle, const FP3& E, const FP3& B, FP& time, double timeStep)
        {
            int thread_id;
#ifdef __USE_OMP__
            thread_id = omp_get_thread_num();
#else
            thread_id = 0;
#endif
            while (time < timeStep)
            {
                FP3 v = particle.getVelocity();
                FP H_eff = sqr(E + (1 / Constants<FP>::lightVelocity()) * VP(v, B))
                    - sqr(SP(E, v) / Constants<FP>::lightVelocity());
                if (H_eff < 0) H_eff = 0;
                H_eff = sqrt(H_eff);
                FP gamma = particle.getGamma();
                FP chi = gamma * H_eff / SchwingerField;
                FP Pdelta = 0.0, dt = 2*timeStep;
                if (chi > 0.0 && coeffPhoton_probability != 0.0)
                {
                    Pdelta = countIntegralParticle(chi);
                    dt = getDtParticle(particle, Pdelta, chi);
                }

                if (dt + time > timeStep)
                {
                    Boris(particle, E, B, timeStep - time);
                    time = timeStep;
                }
                else
                {
                    Boris(particle, E, B, dt);
                    time += dt;
                    FP3 v = particle.getVelocity();
                    FP H_eff = sqr(E + (1 / Constants<FP>::lightVelocity()) * VP(v, B))
                        - sqr(SP(E, v) / Constants<FP>::lightVelocity());
                    if (H_eff < 0) H_eff = 0;
                    H_eff = sqrt(H_eff);
                    FP gamma = particle.getGamma();
                    FP chi_new = gamma * H_eff / SchwingerField;

                    FP delta = Photon_Generator((chi + chi_new)/(FP)2.0);
                    
                    Particle3d NewParticle;
                    NewParticle.setType(Photon);
                    NewParticle.setWeight(particle.getWeight());
                    NewParticle.setPosition(particle.getPosition());
                    NewParticle.setMomentum(delta * particle.getMomentum());

                    this->AvalanchePhotons[thread_id].push_back(NewParticle);
                    this->timeAvalanchePhotons[thread_id].push_back(time);
                    particle.setMomentum((1 - delta) * particle.getMomentum());
                }
            }
        }

        void onePhotonStep(Particle3d& particle, const FP3& E, const FP3& B, FP& time, double timeStep)
        {
            int thread_id;
#ifdef __USE_OMP__
            thread_id = omp_get_thread_num();
#else
            thread_id = 0;
#endif

            FP3 k_ = particle.getVelocity();
            k_ = (1 / k_.norm()) * k_; // normalized wave vector
            FP H_eff = sqrt(sqr(E + VP(k_, B)) - sqr(SP(E, k_)));
            FP gamma = particle.getMomentum().norm()
                / (Constants<FP>::electronMass() * Constants<FP>::lightVelocity());
            FP chi = gamma * H_eff / SchwingerField;

            FP Pdelta = 0.0, dt = 2 * timeStep;
            if (chi > 0.0 && coeffPair_probability != 0.0)
            {
                Pdelta = countIntegralPhoton(chi);
                dt = getDtPhoton(particle, Pdelta, chi, gamma);

            }

            if (dt + time > timeStep)
            {
                particle.setPosition(particle.getPosition()
                    + (timeStep - time) * Constants<FP>::lightVelocity() * k_);
                time = timeStep;
            }
            else
            {
                particle.setPosition(particle.getPosition()
                    + dt * Constants<FP>::lightVelocity() * k_);
                time += dt;
                FP delta = Pair_Generator(chi);

                Particle3d NewParticle;
                NewParticle.setType(Electron);
                NewParticle.setWeight(particle.getWeight());
                NewParticle.setPosition(particle.getPosition());
                NewParticle.setMomentum(delta * particle.getMomentum());

                this->AvalancheParticles[thread_id].push_back(NewParticle);
                this->timeAvalancheParticles[thread_id].push_back(time);

                NewParticle.setType(Positron);
                NewParticle.setMomentum((1 - delta) * particle.getMomentum());
                this->AvalancheParticles[thread_id].push_back(NewParticle);
                this->timeAvalancheParticles[thread_id].push_back(time);

                particle.setP((FP)0.0 * particle.getP());
                time = timeStep;
            }
        }
        
        FP countIntegralParticle(const FP& chi)
        {
            FP a = ((FP)2.0) / ((FP)3.0 * chi);

            FP integral1 = 0;
            FP integral2 = 0;


            if (a < 1e-12)
            {
                integral1 = (FP)5.1984172771983114933 - (FP)18137.993642323330899 * a;
                integral1 *= pow(a, (FP)(1.0) / (FP)(3.0));
            }
            else if (a < 1e-3)
            {
                FP x = pow(a, (FP)(1.0) / (FP)(6.0));
                integral1 = 5.19841743069464487645387850797 +
                    (-0.000018053825036749122040219159830 +
                        (0.00062816225703992768972264921843 +
                            (-0.0093223994163490687907639882599 +
                                (-1.74656176536974939211816042685 +
                                    (-0.235594219314659789134082778933 +
                                        0.327061028774213416726882644466 * x) * x) * x) * x) * x) * x;
                integral1 *= x * x;
            }
            else if (a < 1.0)
            {
                FP x = pow(a, ((FP)(1.0)) / ((FP)(3.0)));
                integral1 = 1.76287424630485660660241735 * 1e-6 +
                    (-0.0000598801738106624577727670992 +
                        (0.000802743773215062663943689639 +
                            (4.38154789085511985519528885788 +
                                (-0.281572569349380384310336886937 +
                                    (0.895653891863824341959839156308 -
                                        0.134897704841544559226394769484 * x) * x) * x) * x) * x) * x;
                integral1 /= ((0.843872941653236121174825417487 + 
                    (-0.057054798244106575113190716375 + 
                        0.464414111160427232538525019242 * x) * x) * x * x);
            }
            else if (a < 4)
            {
                integral1 = 0.10236431861845825264 +
                    (4.4573418932367889106 +
                        (-0.88493986216218777850 +
                            (0.25960997822897173218 +
                                (-0.055661795899496896051 +
                                    (0.0070390038425495395473 -
                                        0.00038680282758132007734 * a) * a) * a) * a) * a) * a;
                integral1 *= pow(a, (FP)(-2.0) / (FP)(3.0));
            }
            else if (a < 100)
            {
                FP x = ((FP)1.0) / a;
                integral1 = 5.23597780390351452 +
                    (-3.222744331519547 +
                        (6.29742855861737 +
                            (-17.43524660816332228 +
                                (46.89295207059157589 +
                                    (-87.51709832603479449758338 +
                                        76.062243906146783344699789 * x) * x) * x) * x) * x) * x;
            }
            else
            {
                FP x = ((FP)1.0) / a;
                integral1 = 5.23598775598298873077107230548 +
                    (-3.22453220308305395661169468026 +
                        (6.39954059064587511538686615114 +
                            (-20.0637559302945579522505446772 +
                                84.6161478096510154145596746648 * x) * x) * x) * x;
            }

            if (a < 1e-3)
            {
                integral2 =
                    0.8664028795330519155244861384644978176518 * pow(a, 1.0 / 3.0) +
                    10.20338035148872944389720486822963776973 * pow(a, 5.0 / 3.0) +
                    11.37153779387130639125888056734653385668 * pow(a, 7.0 / 3.0) +
                    (-2.720699046351326775891117386463233598426 +
                        (-16.96460032938488348769827426970931557466 -
                            4.591179640717863934316260589656706697345 * a) * a) * a;
            }
            else if (a < 4)
            {
                FP t = a * a;
                FP p1 = 10.721735634227961295419410039301150744537629955577 +
                    (3.1068665758080696713212133277483806517693819872828 +
                        (0.22518255645936812947001936658842278197580154122505 +
                            (0.0070981457908672691139275037055531228822662213459040 +
                                (0.00012355857574565104644847829943519078732802492394226 +
                                    (1.3569372024140595131801479874574618752489918107418 * 1e-6 +
                                        (1.0233742472524539571934374951354929985747218340552 * 1e-8 +
                                            (5.6192589175110679089401209364346495908307449422677 * 1e-11 +
                                                (2.3517547519624772649641783723954402272924299437778 * 1e-13 +
                                                    (7.4732788003295848028752816444661730928555234163192 * 1e-16 +
                                                        2.5199767877292697793935365315523079903241129533065 * 1e-18 * 
                                                        t) * t) * t) * t) * t) * t) * t) * t) * t) * t;
                FP p2 = 4.3619451002615210136142189247291646114412066333006 +
                    (0.53806888572279765512786506093776285799130692595714 +
                        (0.022745639260360196230315708847280503764856347206174 +
                            (0.00048478609008443155332306875434206353417875323183142 +
                                (6.2139711935440914661556992333483842213890389804730 * 1e-6 +
                                    (5.3107767139267538162000460941004374413994157565984 * 1e-8 +
                                        (3.2397291785133845863268098122725638677369366881559 * 1e-10 +
                                            (1.4805621798608684459337403067364004591294336758225 * 1e-12 +
                                                (5.2655068576788734519740328960953039075558117454560 * 1e-15 +
                                                    (1.4600683978731546770240400210793190420407004318530 * 1e-17 +
                                                        4.1754686427765817247769862842649983779472220214173 * 1e-20 *
                                                        t) * t) * t) * t) * t) * t) * t) * t) * t) * t;
                FP p3 = 0.86640287953307167358779284429002513122519378768689 +
                    (0.64980215964964131967344556368394889102935543466149 +
                        (0.060918952467505886923142005720464818886544661373976 +
                            (0.0021756768734829432318204947365391972277524448681161 +
                                (0.000040793941570733490593867125526322162594529603228197 +
                                    (4.7069926438452385673387888997826755924774353661428 * 1e-7 +
                                        (3.6773501738593691810668399168843317176831490308209 * 1e-9 +
                                            (2.0735397129198640666107326434756429592099605764322 * 1e-11 +
                                                (8.8485412463143374615380701093381344142964880044152 * 1e-14 +
                                                    (2.8875742154561578697772378331384941634400080057254 * 1e-16 +
                                                        9.4675785148814734522562151731051124873968659647954 * 1e-19 *
                                                        t) * t) * t) * t) * t) * t) * t) * t) * t) * t;
                FP p4 = 10.203380351489012780769226684424763649288872684598 +
                    (1.5305070527224735769012333456788005959584426254853 +
                        (0.071742518097326724374371461406997338039290145372090 +
                            (0.0016305117743927552505020598618101016460351661010823 +
                                (0.000021837211448962242102840630725723229078893264117555 +
                                    (1.9268123765474564257054860604719482642636422353581 * 1e-7 +
                                        (1.2042634137219118467927821316859126535700458205117 * 1e-9 +
                                            (5.6093809853003770270160935319299731278604994152131 * 1e-12 +
                                                (2.0260360527545079572371869231817206147589736101765 * 1e-14 +
                                                    (5.6845234657693501907990135348705608416685575778056 * 1e-17 +
                                                        1.6498017419484357853924900415229540827464506758200 * 1e-19 *
                                                        t) * t) * t) * t) * t) * t) * t) * t) * t) * t;
                FP p5 = 2.7206990463650773295712677753583685969220765603233 +
                    (4.5911796406432008687341397093235187422528594779971 +
                        (0.73786815666846641035074198700499927863037616733811 +
                            (0.038738078122568025201581194023545935326961005805494 +
                                (0.00097957806776202957638023867169463549450918474771544 +
                                    (0.000014431272861671529667921099323302160749543609593385 +
                                        (1.3860786541777686176872491080466454626184663020788 * 1e-7 +
                                            (9.3446897995430291224140740636713642103622728752472 * 1e-10 +
                                                (4.6762037227036676238574617239487881604478613579624 * 1e-12 +
                                                    (1.7330249992940151516415790143624316491177242266850 * 1e-14 +
                                                        6.8211782740989340896898532592041325081051457225634 * 1e-17 *
                                                        t) * t) * t) * t) * t) * t) * t) * t) * t) * t;
                FP p6 = 16.964600329401102757951309576039664796911579641940 +
                    (6.6095845438671531049023635922068737839359517139465 +
                        (0.56525406286108523515168526043051200377148426532569 +
                            (0.019956665146576028807883976050085467096077183658591 +
                                (0.00037848848882886684487878528551834723392826644371587 +
                                    (4.4524300424633283397944801461973636050728790096760 * 1e-6 +
                                        (3.5559314441092138107775713751143602900777989007343 * 1e-8 +
                                            (2.0506941663847479153150985536580737179787468533656 * 1e-10 +
                                                (8.9594392891983095529238980286882022118110001432424 * 1e-13 +
                                                    (2.9549108833342333148746157731369275170435276611438 * 1e-15 +
                                                        1.0345422690524159229324021252792660759968618254786 * 1e-17 *
                                                        t) * t) * t) * t) * t) * t) * t) * t) * t) * t;
                integral2 = pow(t, 7.0 / 6.0) * p1 + pow(t, 11.0 / 6.0) * p2 +
                    pow(t, 1.0 / 6.0) * p3 + pow(t, 5.0 / 6.0) * p4 -
                    pow(t, 1.0 / 2.0) * p5 - t * p6;
            }
            else if (a < 100)
            {
                FP t = pow(a, (FP)(1.0) / (FP)(4.0));

                integral2 = -18.42850180614264514417251180747169259555573144937644239784028126370902 +
                    (93.21237721314488222777689616698495403941636516295620767342606172802269 +
                        (-203.7118110154600987180497321281057176487640593813042577133622177342978 +
                            (252.5174945660205941386056832532073623981519139905190207971273565013861 +
                                (-197.6445837977787615667630787289590732847139300382738282325630359770856 +
                                    (103.1500763595946010601159462623583685481135665922912805189748535358899 +
                                        (-36.64678349946986194393127890576942062898421655029713272175320459118162 +
                                            (8.799130369106253249076661806779624173353280898766252828124392325121153 +
                                                (-1.371360768355051190064627114522550066169432592203585172087251529088909 +
                                                    (0.1255670842498792440178801354356166148844704274545694079188843731314128 -
                                                        0.005137414535995852749595530373871153834004497848623966571298956637843652 *
                                                        t) * t) * t) * t) * t) * t) * t) * t) * t) * t;
                integral2 /= (0.572957795130823208767981548141 * (a * a - a));
            }
            else
            {
                FP x = ((FP)1.0) / a;
                integral2 = (1.745329251994329576923690768488612713442871888541725456097191440171009 +
                    (-12.89812881233221582644677872101088520735335998643726631516362620720443 +
                        (89.59356826904225161541612611574878595673409027847524007965582726211179 +
                            (-668.7918643431519317416848225709347885294334807782286237492250625957854 +
                                (5500.049607627316001946378853216800471232842764317507793778871618035194 -
                                    49936.45920428867757004580008529646421019769989810773723994213800715197 * x) * x) * x) * x) * x) * x * x;
            }

            return integral1 + integral2;
        }

        FP getDtParticle(Particle3d& particle, FP Pdelta, FP chi)
        {
            FP r = -log(random_number_omp());
            r *= ((FP)2.0 * Constants<FP>::pi() * particle.getGamma()) / (sqrt((FP)3) * chi);
            return r / (preFactor*Pdelta);
        }
        
        FP Photon_Generator(FP chi)
        {
            double* g = g_emis;
            FP r = random_number_omp();
            int N = 192;
            FP a = ((FP)2.0) / ((FP)3.0 * chi);
            FP x_target = a;

            FP logA = log(a) / log((FP)1.5);

            int index_a = std::max(std::min((int)std::floor(logA), 29), -30);

            FP delta;
            int index_f = 0;
            FP f1, f2;
            if (r > 1e-1)
            {
                if (r < 0.93)
                {
                    index_f = std::floor((r - 0.1) / 0.83 * 95);
                    f1 = 0.1 + 0.83 / 95.0 * index_f;
                    f2 = f1 + 0.83 / 95.0;
                }
                else if (r < 1.0)
                {
                    index_f = std::min((int)std::floor(-log2(1.0 - (r - 0.93) / 0.07) / 0.2 + 95), N - 1);
                    f1 = (0.07 * (1.0 - pow(2.0, -(index_f - 95) * 0.2)) + 0.93);
                    f2 = (0.07 * (1.0 - pow(2.0, -(index_f - 94) * 0.2)) + 0.93);
                }
                else
                {
                    return 1.0;
                }
            }


            index_f += (index_a + 30) * 192;

            FP x1 = pow(1.5, -index_a);
            FP x2 = pow(1.5, -(index_a + 1));
            x_target = 1.0 / x_target;

            if (r > 1e-1)
            {
                FP z1 = g[index_f] + (g[index_f + 1] - g[index_f]) * (r - f1) / (f2 - f1);
                index_f += N;
                FP z2 = g[index_f] + (g[index_f + 1] - g[index_f]) * (r - f1) / (f2 - f1);

                delta = z1 + (z2 - z1) * (x_target - x1) / (x2 - x1);
            }
            else
            {
                FP z1 = (x2 - x_target) / (x2 - x1) * g[index_f]
                    + (x_target - x1) / (x2 - x1) * g[index_f + N];
                delta = pow(r, 3.0) * z1 / pow(1e-1, 3.0);
            }

            return delta;
        }

        FP countIntegralPhoton(const FP& chi)
        {
            FP a = ((FP)2.0) / ((FP)3.0 * chi);

            FP integral1 = 0;
            FP integral2 = 0;

            FP b = 4 * a;

            if (b < 1e-6)
            {
                FP x = pow(b, 1.0 / 3.0);
                integral1 = 1.5888757318580334808 + (4.621242391684 * 1e-9 +
                    (-1.8138014771310675045 + 0.000258334152383
                        * x) * x) * x;
                integral1 *= x;
            }
            else if (b < 1e-1)
            {
                FP x = pow(b, 1.0 / 3.0);
                integral1 = 1.58887575649 + (-2.88997830078 * 1e-6 +
                    (-1.81369737480 + (-0.00159621376721 +
                        (0.0131013808926 + (0.881783275752 +
                            (-0.100442434842 - 0.224972034631
                                * x) * x) * x) * x) * x) * x) * x;
                integral1 *= x;
            }
            else if (b < 1)
            {
                FP x = pow(b, 1.0 / 3.0);
                integral1 = 0.00222274692761579 + (1.55966451731044 +
                    (0.169055775710587 + (-2.37812332098245 +
                        (2.77990824833067 + (-1.62549030970 +
                            (0.507039151066 - 0.0676910228287
                                * x) * x) * x) * x) * x) * x) * x;
                integral1 *= exp(-b);
            }
            else if (b < 4)
            {

                FP x = pow(b, 1.0 / 3.0);
                integral1 = -0.0863337964364 + (2.12727019609 + 
                    (-1.40279492272 + (0.0617954715806 + 
                        (0.485497729909 + (-0.317700617716 + 
                            (0.0885304524273 - 0.00967875703875 
                                * x) * x) * x) * x) * x) * x) * x;
                integral1 *= exp(-b);
            }
            else if (b < 200)
            {
                FP x = 1.0 / b;
                integral1 = 1.11072184643 + (-0.262362340096 +
                    (0.143672779384 + (0.160152626879 + 
                        (-1.40283063884 + (4.65809948572 +
                            (-8.51231822514 + 6.75456138684
                                * x) * x) * x) * x) * x) * x) * x;
                integral1 *= exp(-b);
            }
            else
                integral1 = 1.11072184643 * exp(-b);
            

            if (a < 0.015625) //2^-6
            {
                FP a1_3 = pow(a, 1.0 / 3.0);
                FP inv_a = 1.0 / a;
                FP a2 = a * a;
                FP a2_3 = a1_3 * a1_3;
                integral2 = 26.4383381757 + (-9.18235928145 +
                    (22.6194671058 + (-5.44139809271 +
                        (2.20691013218 * a1_3
                            - 6.70820546145 * a * a2_3
                            - 8.27591299568 * a2 * a1_3
                            - 9.15670045488 * a2 * a * a2_3)
                        * inv_a) * inv_a) * inv_a) * inv_a;
                integral2 *= a2 * a2;
            }
            else if (a < 8)
            {
                FP integral2_1, integral2_2;
                if (a < 0.5) //integral2_1
                {
                    FP x = pow(a, 1.0 / 3.0);
                    integral2_1 = 1.57027828224594 + (-1.08877277807399 +
                        (-0.179460809334361 + (4.80324684519270 +
                            (-9.74348350227989 + (9.85936580010627 +
                                (-5.290476246 + 1.20361312
                                    * x) * x) * x) * x) * x) * x) * x;
                    integral2_1 *= exp(-4 * a);
                }
                else
                {
                    FP x = pow(a, 1.0 / 6.0);
                    integral2_1 = -1.20066778173054 + (10.2368863985165 +
                        (-25.1115805339200 + (35.9639393894849 +
                            (-29.9577636462431 + (14.7275473913076 +
                                (-3.98478567702 + 0.459485809949
                                    * x) * x) * x) * x) * x) * x) * x;
                    integral2_1 *= exp(-4 * a) / x;
                }

                if (a < 0.125) //integral2_2
                {
                    FP x = pow(a, -1.0 / 6.0);
                    integral2_2 = 4.55073745530463 + (-18.9753607092870 +
                        (32.3758110391415 + (-29.1915770298345 +
                            (15.3907251618537 + (-4.81171065623 +
                                (0.832012023685 - 0.0616171211135
                                    * x) * x) * x) * x) * x) * x) * x;
                    integral2_2 *= exp(-4 * a) * pow(a, -1.0 / 3.0);
                }
                else if (a < 1)
                {
                    FP x = pow(a, -1.0 / 4.0);
                    integral2_2 = -0.311884654115626 + (1.80139927643145 +
                        (-4.30266204865178 + (5.37732685991004 +
                            (-3.50874100787307 + (1.282451108030 +
                                (-0.2517103475662 + 0.02079270065
                                    * x) * x) * x) * x) * x) * x) * x;
                    integral2_2 *= exp(-4 * a);
                }
                else
                {
                    FP x = pow(a, 1.0 / 6.0);
                    integral2_2 = 2.71429586590 + (-16.1475601299 +
                        (39.5087568627 + (-51.5369108153 +
                            (39.5252842521 + (-17.9992925987 +
                                (4.52979557008 - 0.487397088117
                                    * x) * x) * x) * x) * x) * x) * x;
                    integral2_2 *= exp(-4 * a) / a;
                }

                integral2 = integral2_1 - integral2_2;
            }
            else if (a < 50)
            {
                FP x = 1 / a;
                integral2 = 1.11072073452 + (-0.111843403492 +
                    (0.0422690451294 + (-0.0268525026114 +
                        (0.023559419306 + (-0.0241254152815 +
                            (0.0187749900469 + 0.00126008980415
                                * x) * x) * x) * x) * x) * x) * x;
                integral2 *= exp(-4 * a);
            }
            else
                integral2 = 1.11072073452 * exp(-4 * a);


            return std::max(-integral1 * 0.25 + integral2, 0.0);
        }

        FP getDtPhoton(Particle3d& particle, FP Pdelta, FP chi, FP gamma)
        {
            FP r = -log(random_number_omp());
            r *= ((FP)2.0 * Constants<FP>::pi() * gamma) / (sqrt((FP)3) * chi);
            return r / (preFactor * Pdelta);
        }

        FP Pair_Generator(FP chi)
        {
            FP r = random_number_omp();
            if (r < 0.5)
                return 1.0 - Pair_Generator_half(1.0 - r, chi);
            else
                return Pair_Generator_half(r, chi);
        }

        FP Pair_Generator_half(FP r, FP chi)
        {
            double* g = g_pair;
            int N = 128;
            FP a = ((FP)2.0) / ((FP)3.0 * chi);
            FP x_target = a;

            FP logA = log(a) / log((FP)1.5);

            int index_a = std::max(std::min((int)std::floor(logA), 29), -30);

            FP delta;
            int index_f = 0;
            FP f1, f2;
            if (r < 0.95)
            {
                index_f = std::floor((r - 0.5) / 0.45 * 64.0);
                f1 = 0.5 + 0.45 / 64.0 * index_f;
                f2 = f1 + 0.45 / 64.0;
            }
            else if (r < 1.0)
            {
                index_f = std::min((int)std::floor(-log2(1.0 - (r - 0.95) / 0.05) / 0.2 + 64), N - 1);
                f1 = (0.05 * (1.0 - pow(2.0, -(index_f - 64) * 0.2)) + 0.95);
                f2 = (0.05 * (1.0 - pow(2.0, -(index_f - 63) * 0.2)) + 0.95);
            }
            else
            {
                return 1.0;
            }
            

            index_f += (index_a + 30) * N;

            FP x1 = pow(1.5, -index_a);
            FP x2 = pow(1.5, -(index_a + 1));
            x_target = 1.0 / x_target;

            FP z1 = g[index_f] + (g[index_f + 1] - g[index_f]) * (r - f1) / (f2 - f1);
            index_f += N;
            FP z2 = g[index_f] + (g[index_f + 1] - g[index_f]) * (r - f1) / (f2 - f1);

            delta = z1 + (z2 - z1) * (x_target - x1) / (x2 - x1);
        
            return delta;
        }

        void operator()(ParticleProxy3d* particle, ValueField field, FP timeStep)
        {}

        void operator()(Particle3d* particle, ValueField field, FP timeStep)
        {
            ParticleProxy3d particleProxy(*particle);
            this->operator()(&particleProxy, field, timeStep);
        }
    private:

        FP random_number_omp()
        {
            FP rand_n;
#pragma omp critical
            rand_n = distribution(rand_generator);
            return rand_n;
        }

        
        std::default_random_engine rand_generator;
        std::uniform_real_distribution<FP> distribution;

        vector<vector<FP>> timeAvalanchePhotons, timeAvalancheParticles;
        vector<vector<Particle3d>> AvalanchePhotons, AvalancheParticles;
        vector<vector<Particle3d>> afterAvalanchePhotons, afterAvalancheParticles;

        inline static const uint64_t int_g_emis[] = {
#include "QED_emis.in" 
        };

        inline static const uint64_t int_g_pair[] = {
#include "QED_pair_half.in" 
        };

        FP SchwingerField;
        FP preFactor;
        FP coeffPhoton_probability, coeffPair_probability;

        double * g_emis, * g_pair;
    };

    typedef Scalar_Fast_QED<YeeGrid> Scalar_Fast_QED_Yee;
    typedef Scalar_Fast_QED<PSTDGrid> Scalar_Fast_QED_PSTD;
    typedef Scalar_Fast_QED<PSATDGrid> Scalar_Fast_QED_PSATD;
    typedef Scalar_Fast_QED<AnalyticalField> Scalar_Fast_QED_Analytical;
}