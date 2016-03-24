library(plotly)

py <- plotly(username="r_user_guide", key="mw5isa4yqp")  # open plotly connection

pp <- function (n,r=4) {
  x <- seq(-r*pi, r*pi, len=n)
  df <- expand.grid(x=x, y=x)
  df$r <- sqrt(df$x^2 + df$y^2)
  df$z <- cos(df$r^2)*exp(-df$r/6)
  df
}
p <- ggplot(pp(20), aes(x=x,y=y))

p <- p + geom_tile(aes(fill=z))

py$ggplotly(p)
