library(rethinking)
library(tidyverse)
library(rstan)
library(gridExtra)
source("functions/functions.R")

# Model
m <- readRDS(file="stanout/mogK3.rda")

param <- c("theta_simple", 
           "theta_complex" 
           )

#round(summary(m, param)$summary,2)
#traceplot(m, param)

pdiffs <- as.data.frame(m, pars = param) %>%
   transmute(
         p1 = `theta_complex[1]` - `theta_simple[1]`,
         p2 = `theta_complex[2]` - `theta_simple[2]`,
         p3 = `theta_complex[3]` - `theta_simple[3]` 
         ) %>%
    gather(Parameter, value) %>% as.tibble()

pdiffs %>%
    group_by(Parameter) %>%
    dplyr::summarise( 
      mu = dmode(value),
      Lower = HPDI(value, prob = .95)[1],
      Upper = HPDI(value, prob = .95)[2]
      )


pdiffs$Parameter <- as.factor(pdiffs$Parameter)
levels(pdiffs$Parameter) <- c("K_1","K_2","K_3")

vnames <-list(
  'p1' = quote(italic("K")[1]),
  'p2' = quote(italic("K")[2]),
  'p3' = quote(italic("K")[3]))

vlabeller <- function(variable,value){
  return(vnames[value])
}

# BPS slides

p <- pdiffs %>% 
  ggplot(aes(x=value, fill = Parameter)) +
  geom_histogram(color="white", bins = 80, alpha =.35) +
  scale_fill_manual(
    values = c( "darkred", "darkblue", "forestgreen"),
    breaks = c("p1", "p2", "p3")) +
  facet_grid(Parameter~., 
             scales = "free_y", 
             space = "free",
             labeller = vlabeller
  ) +
  theme_minimal() +
  xlab(bquote(Delta*hat(theta)~"(conjoined NP"-"simple NP)")) + #~(conjoined~NP~-~simple~NP)
  #  ylab("Posterior probability") +
#  ylab("") +
  ylab("Posterior probability") +
#  ggtitle("Mixing proportion") +
  geom_vline(aes(xintercept=0), linetype = "dotted", color = "grey40") +
  scale_x_continuous(limits = c(-.37,.37), breaks = seq(-.3,.3,.15)) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
 #       plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size =10),
        panel.grid.minor = element_blank(),
        text=element_text(family="ArialMT"));p 

#cairo_pdf("slides-ntu-bps/gfx/results_p.pdf", width = 6.5, height = 5,  bg = "transparent")
cairo_pdf("slides-tilburg-coll/gfx/results_p.pdf", width = 6.5, height = 5,  bg = "transparent")
p
dev.off()

param <- c("theta_simple", 
           "theta_complex" ,
           "NP0", "NP",  "NP2"  
#           "sigmap_e[2]", # K1
 #           "sigma_e" , # K2
  #         "sigmap_e[1]" # K3
)


ps <- as.data.frame(m, pars = param) %>%
  transmute(
    #          sigma0 = `sigmap_e[2]`,
#         sigma1 = exp(sigma_e),
#         sigma2 = `sigmap_e[1]`,
#         sigma0 = exp(sigma0),
#         sigma2 = exp(sigma2),
         p_simple_1 = `theta_simple[1]`,
         p_simple_2 = `theta_simple[2]`,
         p_simple_3 = `theta_simple[3]`,
         p_conjoined_1 = `theta_complex[1]`,
         p_conjoined_2 = `theta_complex[2]`,
         p_conjoined_3 = `theta_complex[3]`
  ) %>%
  gather(Parameter, value) %>% 
  as.tibble(); ps

ps$Parameter <- as.factor(ps$Parameter)

ps$Parameter <- factor(ps$Parameter, levels = c(  "p_conjoined_1",  "p_conjoined_2", "p_conjoined_3", "p_simple_1",    "p_simple_2",   "p_simple_3"   )) 

vnames <-list(
#  'NP0' = quote(hat(mu)[1]),
#  'NP1' = quote(hat(mu)[2]),
#  'NP2' = quote(hat(mu)[3]),
#  "sigma0" = quote(hat(sigma)[1]),
#  "sigma1" = quote(hat(sigma)[2]),
#  "sigma2" = quote(hat(sigma)[3]),
  "p_conjoined_1" = quote(hat(theta)["conjoined"[1]]),
  "p_conjoined_2" = quote(hat(theta)["conjoined"[2]]),
  "p_conjoined_3" = quote(hat(theta)["conjoined"[3]]),
  "p_simple_1" = quote(hat(theta)["simple"[1]]),
  "p_simple_2" = quote(hat(theta)["simple"[2]]),
  "p_simple_3" = quote(hat(theta)["simple"[3]])
)


vlabeller <- function(variable,value){
  return(vnames[value])
}

p <- ps %>% 
  ggplot(aes(x=value)) +
  geom_histogram(color="white", bins = 20, alpha =.35, fill = "forestgreen") +
  facet_wrap(~Parameter, 
             scales = "free",
             nrow = 2,
             labeller = vlabeller
  ) +
  theme_minimal() +
  xlab("") + #~(conjoined~NP~-~simple~NP)
  ylab("") +
  ggtitle(bquote("Mixing proportion"~hat(theta)),
          subtitle = "by NP and mixture component") +
  theme(strip.text = element_text(angle = 0, size = 12),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10),
 #       plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 9),
        panel.grid.minor = element_blank());p

cairo_pdf("~/Documents/poster_cuny2019/gfx/mixp.pdf", width = 4.5, height = 4.45,  bg = "transparent")
p
dev.off()



# CUNY abstract:
p <- pdiffs %>% 
  filter(Parameter == "K_3") %>%
  ggplot(aes(x=value, fill = Parameter)) +
  geom_histogram(color="white", bins = 20, alpha =.35) +
  scale_fill_manual(
  values = c("forestgreen"),
  breaks = c("p3")) +
  theme_minimal() +
  ggtitle(bquote("Difference"~Delta*theta[3]), subtitle = "between conjoined and simple NP")+
  xlab(bquote(Delta*theta[3])) + #~(conjoined~NP~-~simple~NP)
  ylab(" ") +
##  geom_vline(aes(xintercept=0), linetype = "dotted", color = "grey40") +
#  scale_x_continuous(limits = c(-.37,.37), breaks = seq(-.3,.3,.15)) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 8),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size = 9),
        axis.title = element_text(size = 10),
        panel.grid.minor = element_blank());p 
  
cairo_pdf("~/Documents/poster_cuny2019/gfx/difmixp.pdf", width = 3, height = 4.3,  bg = "transparent")
p
dev.off()
