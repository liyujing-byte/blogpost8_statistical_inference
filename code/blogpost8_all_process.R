# Blog Post 8: Simulation and Statistical Inference
# Run from the project root: source("code/blogpost8_all_process.R")
# Base R only; all data are simulated, not real household observations.
set.seed(640008)
for (d in c("code", "data/raw", "data/processed", "data/results", "figures/post", "figures/evidence")) dir.create(d, showWarnings = FALSE, recursive = TRUE)
mu <- 100; sigma <- 20; n <- 80; B <- 2000; M <- 1000
x <- rnorm(n, mu, sigma)
write.csv(data.frame(spending = x), "data/raw/sample.csv", row.names = FALSE)
boot_mean <- function(z, b = B) replicate(b, mean(sample(z, replace = TRUE)))
bm <- boot_mean(x)
ci <- quantile(bm, c(.025, .975))
write.csv(data.frame(bootstrap_mean = bm), "data/results/bootstrap_draws.csv", row.names = FALSE)
# Each coverage trial draws a NEW sample, then bootstraps that sample.
intervals <- t(replicate(M, {
  z <- rnorm(n, mu, sigma)
  quantile(boot_mean(z), c(.025, .975))
}))
hit <- intervals[,1] <= mu & intervals[,2] >= mu
coverage <- mean(hit)
mcse <- sqrt(coverage * (1 - coverage) / M)
write.csv(data.frame(lower=intervals[,1], upper=intervals[,2], covers=hit),
          "data/results/coverage.csv", row.names = FALSE)
# Two-sided test with KNOWN population SD: H0: mean = 95.
# Simulate Z under H0; do not bootstrap the uncentered observed sample.
mu0 <- 95; S <- 50000; alpha <- .05
zobs <- (mean(x) - mu0) / (sigma / sqrt(n))
znull <- replicate(S, (mean(rnorm(n, mu0, sigma)) - mu0) / (sigma / sqrt(n)))
critical <- unname(quantile(abs(znull), 1-alpha))
p_mc <- (1 + sum(abs(znull) >= abs(zobs))) / (S + 1)
p_exact <- 2 * pnorm(-abs(zobs))
write.csv(data.frame(z=znull), "data/results/null_draws.csv", row.names=FALSE)
summary_table <- data.frame(sample_mean=mean(x), lower=ci[1], upper=ci[2],
                            coverage=coverage, coverage_mcse=mcse, z_observed=zobs,
                            critical_abs_z=critical, simulated_p=p_mc, exact_p=p_exact)
write.csv(summary_table, "data/results/summary.csv", row.names=FALSE)
# Readable figures, saved automatically.
png("figures/01_bootstrap.png", 1200, 760, res=140)
par(mar=c(5,5,4,2))
hist(bm, breaks=35, col="#8fc9df", border="white",
     main="One sample, many resamples", xlab="Bootstrap mean spending ($)")
abline(v=ci, col="#d45b41", lwd=3, lty=2)
abline(v=mean(x), col="#194d66", lwd=3)
legend("topright", c("Sample mean", "95% percentile endpoints"),
       col=c("#194d66", "#d45b41"), lty=c(1,2), lwd=3, bty="n")
dev.off()
png("figures/02_coverage.png", 1200, 860, res=140)
par(mar=c(5,5,4,2)); k <- 100
plot(range(intervals[1:k,]), c(1,k), type="n", xlab="95% bootstrap interval for mean spending ($)",
     ylab="Independent sample", main=sprintf("First 100 intervals; all %s trials: %.1f%% coverage",M,100*coverage))
segments(intervals[1:k,1],1:k,intervals[1:k,2],1:k,
         col=ifelse(hit[1:k],"#428da8","#d45b41"),lwd=2)
abline(v=mu,lty=2,lwd=2)
legend("topright",c("Contains true mean", "Misses true mean"),
       col=c("#428da8","#d45b41"),lwd=2,bty="n")
dev.off()
png("figures/03_null_test.png",1200,760,res=140)
par(mar=c(5,5,4,2)); lim <- max(4.5,abs(zobs)+.5)
h <- hist(znull,breaks=seq(-8,8,.16),plot=FALSE)
plot(h,freq=FALSE,xlim=c(-lim,lim),col=ifelse(abs(h$mids)>=critical,"#e99781","#8fc9df"),
     border="white", main="What would we observe if H0 were true?", xlab="Z statistic under H0: mean spending = $95")
abline(v=c(-critical,critical),col="#b74732",lty=2,lwd=2)
abline(v=c(-abs(zobs),abs(zobs)),col="#194d66",lwd=3)
legend("topright",c("5% rejection boundaries", "Observed |Z| boundaries"),
       col=c("#b74732","#194d66"),lty=c(2,1),lwd=2,bty="n")
dev.off()
# Narrative incorporates computed results, avoiding manually copied numbers.
paragraphs <- c(
  "# Learning Statistical Inference by Repeating an Experiment",
  "## Question and motivation",
  "How can one small sample tell us something about a population? A sample average changes whenever we collect new observations. Statistical inference quantifies this uncertainty. Here I simulate independent household spending from a normal population with mean $100 and standard deviation $20. These numbers are assumptions for a teaching experiment, not estimates from real households. Knowing the population truth lets us check whether our methods work.",
  "## From one sample to a bootstrap distribution",
  sprintf("I draw %s observations and obtain a sample mean of $%.2f. To bootstrap, I draw %s observations with replacement from this same sample, calculate their mean, and repeat %s times. Replacement matters: an observation can appear several times in one resample. This approximates the sampling variability of the mean using the observed data as an empirical population.",n,mean(x),n,B),
  sprintf("The first figure displays the resulting bootstrap distribution. Its 2.5th and 97.5th percentiles give a 95%% percentile interval of [$%.2f, $%.2f]. The interval describes uncertainty about the population mean, not the range containing 95%% of individual households. A narrow distribution means greater precision; a larger sample generally makes the distribution narrower.",ci[1],ci[2]),
  "![Bootstrap distribution](../01_bootstrap.png)",
  "## Does 95% really mean 95%?",
  sprintf("To investigate coverage, I repeat the entire experiment %s times. Every repetition starts with a fresh sample from the population, constructs its own bootstrap distribution, and records whether its interval includes the fixed true mean of $100. The second figure shows the first 100 intervals; orange intervals miss the vertical truth line. Across all repetitions, coverage is %.1f%%, with a Monte Carlo standard error of %.2f percentage points.",M,100*coverage,100*mcse),
  "The 95% refers to the long-run success rate of the interval procedure. After a particular interval is calculated, both endpoints are fixed, as is the population mean: the interval either contains it or does not. Observed coverage need not equal exactly 95%. Finite samples, bootstrap approximation and simulation randomness all matter. Percentile intervals can under-cover, especially for small or strongly skewed samples; this experiment does not establish a universal guarantee.",
  "![Repeated interval coverage](../02_coverage.png)",
  "## Hypothesis testing: a world where the null is true",
  "Now consider H0: mean spending equals $95, against a two-sided alternative. For this demonstration I treat the population standard deviation of $20 as known. The test statistic is Z = (sample mean - 95)/(20/sqrt(n)). If the standard deviation were unknown, the testing procedure would need to account for estimating it, for example through a t test under normality.",
  sprintf("I simulate %s new samples under H0 to build the reference distribution. The observed statistic is %.3f. The outer 5%% of absolute simulated statistics defines a two-sided rejection region, approximately |Z| > %.3f. The orange histogram tails show this rejection region. Dark lines mark the observed absolute statistic on both sides, separating outcomes at least as extreme as our observation.",S,zobs,critical),
  sprintf("The Monte Carlo p-value is %.4f; the normal-theory value is %.4f. The simulated calculation counts both tails and uses a plus-one correction to avoid reporting zero from a finite simulation. At the 5%% level, we %s H0. A p-value is the probability, assuming H0 and the model, of a statistic at least this extreme. It is not the probability that H0 is true and does not measure economic importance.",p_mc,p_exact,ifelse(p_mc<.05,"reject","do not reject")),
  "![Null distribution and rejection regions](../03_null_test.png)",
  "## Try it and takeaway",
  "Open interactive.html to vary sample size and confidence level. Its seeded normal-population experiment recomputes bootstrap distributions and repeated-sample coverage in your browser. Compare interval width and coverage across settings; small simulation counts add visible noise. The R analysis remains the reproducible source for the reported figures and tables.",
  "Bootstrap intervals quantify estimation uncertainty; null simulations calibrate tests. Both rely on a sampling model and repeated experiments. Independence and representative sampling are substantive assumptions: repeating a biased sampling process cannot repair it.",
  "## Data and reproducibility",
  "All observations are synthetic and generated with base R using seed 640008. Run source('code/blogpost8_all_process.R') from the project root. See README.md for outputs and assumptions. No external dataset or additional R package is required."
)
writeLines(paragraphs,"figures/post/blogpost8.md")
# Self-contained published-style HTML with embedded figures.
body <- vapply(paragraphs,function(p) {
  if(startsWith(p,"![")) {
    path <- sub(".*\\((.*)\\)","\\1",p)
    return(paste0('<img src="',basename(path),'" alt="Simulation result">'))
  }
  if(startsWith(p,"# ")) return(paste0('<h1>',substring(p,3),'</h1>'))
  if(startsWith(p,"## ")) return(paste0('<h2>',substring(p,4),'</h2>'))
  paste0('<p>',p,'</p>')
},character(1))
writeLines(c('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Blog Post 8</title><style>body{max-width:850px;margin:50px auto;padding:0 24px;font:18px/1.7 Georgia;color:#183747}h1,h2{font-family:Arial}img{width:100%}a{color:#176f91}</style>',body,
             '<p><a href="figures/interactive.html">Open the interactive inference lab</a></p></html>'),"figures/index.html")
writeLines(capture.output(sessionInfo()),"data/results/sessionInfo.txt")
print(summary_table)
cat("\nComplete: index.html, blogpost8.md, results/ and data/ generated.\n")

# Standalone interactive companion (generated with the analysis).
writeLines("<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Inference lab</title><style>body{max-width:900px;margin:40px auto;padding:20px;font:17px/1.6 Arial;color:#183747}canvas{width:100%;border:1px solid #ddd}label{margin-right:20px}button{padding:10px;background:#176f91;color:white;border:0}input{width:65px}</style><h1>Bootstrap inference lab</h1><p>Synthetic independent normal spending: mean $100, SD $20. Change settings and rerun. Seed 640008 resets on every run for reproducibility.</p><label>Sample size <input id=\"n\" type=\"number\" min=\"10\" max=\"200\" value=\"80\"></label><label>Confidence % <input id=\"cl\" type=\"number\" min=\"80\" max=\"99\" value=\"95\"></label><button id=\"run\">Run simulation</button><p id=\"out\" aria-live=\"polite\"></p><h2>Bootstrap means from one sample</h2><canvas id=\"hist\" width=\"850\" height=\"300\"></canvas><h2>First 60 repeated-sample intervals</h2><p>Blue contains $100; orange misses. Black line is the true mean.</p><canvas id=\"cov\" width=\"850\" height=\"350\"></canvas><p>500 resamples per interval; 200 independent samples for coverage. Coverage is approximate and has Monte Carlo noise. Confidence describes the procedure, not the probability assigned to one fixed interval.</p><a href=\"index.html\">Back to blog</a><script>\nlet seed;function u(){seed=(1664525*seed+1013904223)>>>0;return (seed+.5)/4294967296}function normal(){return 100+20*Math.sqrt(-2*Math.log(u()))*Math.cos(2*Math.PI*u())}function avg(a){return a.reduce((s,v)=>s+v,0)/a.length}function quant(a,p){let j=(a.length-1)*p,k=Math.floor(j);return a[k]+(a[Math.min(k+1,a.length-1)]-a[k])*(j-k)}function boots(x){let b=[];for(let j=0;j<500;j++){let s=0;for(let i=0;i<x.length;i++)s+=x[Math.floor(u()*x.length)];b.push(s/x.length)}return b.sort((a,b)=>a-b)}\nfunction run(){seed=640008;let n=Math.min(200,Math.max(10,Math.round(+document.getElementById('n').value||80))),cl=Math.min(99,Math.max(80,+document.getElementById('cl').value||95)),a=(1-cl/100)/2;let x=Array.from({length:n},normal),b=boots(x),lo=quant(b,a),hi=quant(b,1-a);let ints=[],hits=0;for(let k=0;k<200;k++){let z=boots(Array.from({length:n},normal)),l=quant(z,a),h=quant(z,1-a);ints.push([l,h]);if(l<=100&&h>=100)hits++}document.getElementById('out').textContent='Sample mean: $'+avg(x).toFixed(2)+'; '+cl+'% bootstrap CI: [$'+lo.toFixed(2)+', $'+hi.toFixed(2)+']; coverage: '+hits+'/200 = '+(hits/2).toFixed(1)+'% (Monte Carlo SE about '+(100*Math.sqrt(hits/200*(1-hits/200)/200)).toFixed(1)+' percentage points).';let c=document.getElementById('hist').getContext('2d');c.clearRect(0,0,850,300);let mn=b[0],mx=b[b.length-1],bins=Array(25).fill(0),scale=v=>50+(v-mn)/(mx-mn)*750;b.forEach(v=>bins[Math.min(24,Math.floor((v-mn)/(mx-mn)*25))]++);let top=Math.max(...bins);bins.forEach((v,i)=>{c.fillStyle='#8fc9df';c.fillRect(50+i*30,250-v/top*200,28,v/top*200)});c.strokeStyle='#c45036';[lo,hi].forEach(v=>{c.beginPath();c.moveTo(scale(v),30);c.lineTo(scale(v),250);c.stroke()});c.fillStyle='#183747';c.fillText('Bootstrap mean ($): '+mn.toFixed(1),50,280);c.fillText(mx.toFixed(1),770,280);c.fillText('Count',10,20);let d=document.getElementById('cov').getContext('2d');d.clearRect(0,0,850,350);let left=Math.min(...ints.slice(0,60).map(v=>v[0]),100)-1,right=Math.max(...ints.slice(0,60).map(v=>v[1]),100)+1,s=v=>50+(v-left)/(right-left)*750;ints.slice(0,60).forEach((v,i)=>{d.strokeStyle=v[0]<=100&&v[1]>=100?'#428da8':'#d45b41';d.beginPath();d.moveTo(s(v[0]),20+5*i);d.lineTo(s(v[1]),20+5*i);d.stroke()});d.strokeStyle='black';d.beginPath();d.moveTo(s(100),10);d.lineTo(s(100),320);d.stroke();d.fillStyle='#183747';d.fillText('Mean spending ($): '+left.toFixed(1),50,340);d.fillText(right.toFixed(1),770,340)}document.getElementById('run').onclick=run;run();</script></html>", "figures/interactive.html")

front <- '\ntitle: "My Eighth Blog Post"\nsubtitle: "Learning Statistical Inference by Repeating an Experiment"\ndate: "2026-11-02"\ncategories: [data science, simulation, statistical inference]\nimage: ../01_bootstrap.png\nresources:\n  - ../interactive.html\n'
writeLines(c("---", front, "---", "", paragraphs[-1], "", "[Complete R code, data, and replication evidence](https://github.com/liyujing-byte/blogpost8_statistical_inference)", "", "```{=html}", '<iframe src="../interactive.html" title="Interactive bootstrap inference lab" style="width:100%;height:1050px;border:1px solid #dde4e9;"></iframe>', "```", "", "[Open the interactive lab](../interactive.html)"), "figures/post/blogpost8.qmd")
write.csv(data.frame(observation=seq_along(x), spending=x), "data/processed/analysis_data.csv", row.names=FALSE)
checks <- data.frame(check=c("sample size","no missing data","bootstrap count","coverage count","null count","ordered CI","valid p-value","coverage near 95 percent"),passed=c(length(x)==n,!anyNA(x),length(bm)==B,nrow(intervals)==M,length(znull)==S,ci[1]<ci[2],p_mc>0 && p_mc<=1,abs(coverage-.95)<.03))
stopifnot(all(checks$passed))
write.csv(checks,"figures/evidence/data_checks.csv",row.names=FALSE)
writeLines("PASS: all eight numerical checks", "figures/evidence/output_validation.log")
paths <- list.files("data", recursive=TRUE, full.names=TRUE)
write.csv(data.frame(path=paths,md5=unname(tools::md5sum(paths))),"figures/evidence/data_result_checksums.csv",row.names=FALSE)
